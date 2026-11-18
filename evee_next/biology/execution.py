"""Explicitly authorized, bounded retrieval and pinned local CPU execution.

Policy limits never grant permission. Transports may be injected for offline
verification; default HTTP handling disables automatic redirects and retries.
"""
from contextlib import contextmanager
from dataclasses import dataclass, replace
from email.policy import HTTP
import hashlib
import json
import os
from pathlib import Path
import re
import selectors
import signal
import subprocess
import threading
import time
from types import MappingProxyType
from urllib.error import HTTPError
from urllib.parse import urljoin, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

from .evidence_store import decision_from_dict
from .routing import EvidenceDecision, canonical_bytes, require_sha


class BudgetExceeded(RuntimeError): pass
class IntegrityError(ValueError): pass
class PendingAuthorization(PermissionError): pass


def require_current_policy(policy,*,fetcher=None,budget=None):
    """An injected executor cannot retain stronger permissions than this call."""
    if fetcher is not None:
        if getattr(fetcher,'policy',policy)!=policy:
            raise PendingAuthorization('Configured fetcher policy differs from current execution policy')
        budget=getattr(fetcher,'budget',budget)
    if budget is not None and getattr(budget,'policy',policy)!=policy:
        raise PendingAuthorization('Configured execution budget differs from current policy')


class ExecutionBudget:
    """One invocation's cumulative accounting, shared by every provider."""
    def __init__(self, policy, *, clock=time.monotonic):
        self.policy, self.clock = policy, clock
        self.started = clock()
        self.bytes_fetched = self.requests = self.redirects = self.retries = 0
        self.variant_bytes = {}
        self.runner_counts = {}
        self._lock, self._runner = threading.RLock(), threading.Lock()
        self._retrieval = threading.Lock()

    def remaining_seconds(self):
        remaining = self.policy.invocation_seconds - (self.clock() - self.started)
        if remaining <= 0: raise BudgetExceeded('Invocation time budget exhausted')
        return remaining

    def remaining_bytes(self, variant_key):
        with self._lock:
            return min(self.policy.invocation_bytes-self.bytes_fetched,
                       self.policy.per_variant_bytes-self.variant_bytes.get(variant_key, 0))

    def consume_bytes(self, size, variant_key):
        if type(size) is not int or size < 0: raise ValueError('Nonnegative integer byte charge required')
        if not isinstance(variant_key, str) or not variant_key: raise ValueError('Variant accounting key required')
        with self._lock:
            self.remaining_seconds()
            if size > self.remaining_bytes(variant_key): raise BudgetExceeded('Fetched byte budget exhausted')
            self.bytes_fetched += size
            self.variant_bytes[variant_key] = self.variant_bytes.get(variant_key, 0) + size

    def retry(self):
        with self._lock:
            self.remaining_seconds()
            if self.retries >= self.policy.transport_retries: return False
            self.retries += 1
            return True

    @contextmanager
    def cpu_runner(self, runner_id):
        self.remaining_seconds()
        if not self._runner.acquire(blocking=False): raise BudgetExceeded('Only one CPU runner may execute at a time')
        try:
            self.runner_counts[runner_id] = self.runner_counts.get(runner_id, 0) + 1
            yield min(self.policy.task_seconds, self.remaining_seconds())
        finally: self._runner.release()


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl): return None


class _HTTPResponse:
    def __init__(self, response):
        self.response = response
        self.status, self.headers = response.code, dict(response.headers.items())
        status_line = f'HTTP/1.1 {response.code} {response.reason}\r\n'.encode()
        self.header_bytes = len(status_line) + len(response.headers.as_bytes(policy=HTTP)) + 2
    def set_timeout(self, remaining):
        # urllib's connected socket remains reachable through the response's
        # buffered reader. Refresh its timeout before every bounded read.
        if self.response.fp is not None:
            self.response.fp.raw._sock.settimeout(remaining)
    def iter_bytes(self, chunk_size):
        while True:
            data = self.response.read(chunk_size)
            if not data: break
            yield data
    def read(self, size): return self.response.read(size)
    def close(self): self.response.close()


def _transport(url, *, headers, timeout):
    try: response = build_opener(_NoRedirect()).open(Request(url, headers=headers), timeout=timeout)
    except HTTPError as exc: response = exc
    return _HTTPResponse(response)


class BoundedFetcher:
    def __init__(self, policy, budget, *, transport=None):
        if budget.policy != policy: raise ValueError('Fetcher and budget policy must match')
        self.policy, self.budget, self.transport = policy, budget, transport or _transport

    def _url(self, url):
        parts = urlsplit(url)
        if (parts.scheme != 'https' or parts.hostname not in self.policy.allowed_hosts
                or parts.username is not None or parts.password is not None
                or parts.port not in (None, 443) or parts.fragment
                or any(ord(c) < 32 for c in url)):
            raise IntegrityError('Only credential-free HTTPS URLs on explicitly allowed hosts are supported')
        return parts

    def fetch(self, url, *, variant_key, expected_etag=None, byte_range=None):
        # Serialize streaming reservations within an invocation. Two readers
        # cannot both spend the same remaining bytes before charging them.
        if not self.budget._retrieval.acquire(timeout=self.budget.remaining_seconds()):
            raise BudgetExceeded('Invocation time budget exhausted waiting for retrieval')
        try:
            return self._fetch(url, variant_key=variant_key, expected_etag=expected_etag, byte_range=byte_range)
        finally:
            self.budget._retrieval.release()

    def _fetch(self, url, *, variant_key, expected_etag=None, byte_range=None):
        if self.policy.mode != 'bounded-retrieval' or not self.policy.authorization_refs:
            raise PendingAuthorization('Current retrieval authorization and bounded policy required')
        if expected_etag is not None and (not re.fullmatch(r'"[^"\r\n]+"', expected_etag)):
            raise IntegrityError('A strong exact ETag pin is required')
        if byte_range is not None:
            if (not isinstance(byte_range, tuple) or len(byte_range) != 2
                    or any(type(v) is not int for v in byte_range)
                    or not 0 <= byte_range[0] <= byte_range[1] or expected_etag is None):
                raise IntegrityError('Inclusive byte range and strong source pin required')
        redirects = 0
        deadline = self.budget.clock() + self.policy.task_seconds
        def remaining_time():
            remaining = min(deadline-self.budget.clock(), self.budget.remaining_seconds())
            if remaining <= 0: raise BudgetExceeded('Retrieval task time budget exhausted')
            return remaining
        while True:
            parts = self._url(url)
            remaining_time()
            headers = {'Accept-Encoding': 'identity', 'User-Agent': 'EVEE-qualified-evidence/1',
                       'Connection': 'close'}
            if expected_etag: headers['If-Match'] = expected_etag
            if byte_range: headers['Range'] = f'bytes={byte_range[0]}-{byte_range[1]}'
            target = (parts.path or '/') + ('?' + parts.query if parts.query else '')
            wire_request = f'GET {target} HTTP/1.1\r\nHost: {parts.netloc}\r\n' + ''.join(f'{k}: {v}\r\n' for k,v in headers.items()) + '\r\n'
            self.budget.consume_bytes(len(wire_request.encode()), variant_key)
            self.budget.requests += 1
            response = None
            try:
                response = self.transport(url, headers=headers, timeout=remaining_time())
                remaining_time()
                normalized = {str(k).lower(): str(v) for k,v in response.headers.items()}
                wire_headers = f'HTTP/1.1 {response.status}\r\n' + ''.join(f'{k}: {v}\r\n' for k,v in response.headers.items()) + '\r\n'
                self.budget.consume_bytes(getattr(response, 'header_bytes', len(wire_headers.encode())), variant_key)
                redirect = response.status in (301,302,303,307,308)
                transient = response.status in (408,429,500,502,503,504)
                if not redirect and not transient:
                    if expected_etag is not None and normalized.get('etag') != expected_etag:
                        raise IntegrityError('Source ETag pin mismatch')
                    if normalized.get('content-encoding', 'identity').lower() != 'identity':
                        raise IntegrityError('Encoded payload cannot satisfy byte integrity')
                    if byte_range:
                        match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', normalized.get('content-range',''))
                        if (response.status != 206 or match is None
                                or tuple(map(int, match.groups()[:2])) != byte_range
                                or int(match.group(3)) <= byte_range[1]):
                            raise IntegrityError('Range retrieval requires exact partial response; whole-file fallback refused')
                    elif response.status != 200:
                        raise IntegrityError(f'Unexpected HTTP status {response.status}')
                body = bytearray()
                declared_size = None
                if 'content-length' in normalized:
                    try: declared_size = int(normalized['content-length'])
                    except ValueError as exc: raise IntegrityError('Invalid content length') from exc
                    if declared_size < 0: raise IntegrityError('Invalid content length')
                if byte_range and not redirect and not transient:
                    range_size = byte_range[1]-byte_range[0]+1
                    if declared_size is not None and declared_size != range_size:
                        raise IntegrityError('Partial response declared length mismatch')
                while True:
                    remaining_time()
                    remaining = self.budget.remaining_bytes(variant_key)
                    if remaining <= 0:
                        if declared_size == len(body): break
                        raise BudgetExceeded('Fetched byte budget exhausted')
                    if hasattr(response, 'set_timeout'): response.set_timeout(remaining_time())
                    chunk = response.read(min(65536, remaining))
                    if not isinstance(chunk, bytes): raise IntegrityError('Transport must stream bytes')
                    if not chunk: break
                    self.budget.consume_bytes(len(chunk), variant_key)
                    body.extend(chunk)
                    remaining_time()
                    if declared_size is not None and len(body) > declared_size:
                        raise IntegrityError('Response exceeds declared content length')
                if declared_size is not None and len(body) != declared_size:
                    raise IntegrityError('Content length mismatch')
                if transient: raise OSError('Transient HTTP response')
                if redirect:
                    self.budget.redirects += 1
                    redirects += 1
                    if redirects > 5 or not normalized.get('location'): raise IntegrityError('Redirect limit or missing location')
                    url = urljoin(url, normalized['location'])
                    continue
                if byte_range and len(body) != byte_range[1]-byte_range[0]+1:
                    raise IntegrityError('Partial response payload length mismatch')
                return bytes(body)
            except (OSError, TimeoutError):
                if not self.budget.retry(): raise
            finally:
                if response is not None: response.close()


@dataclass(frozen=True)
class RunnerSpec:
    runner_id: str
    command: tuple[str, ...]
    file_sha256: dict
    authorization_refs: tuple[str, ...]
    max_output_bytes: int = 1024**2
    device: str = 'cpu'

    def __post_init__(self):
        if not self.runner_id or not self.command or not self.authorization_refs:
            raise ValueError('Qualified runner ID, command and current authorization bindings required')
        if self.device != 'cpu' or type(self.max_output_bytes) is not int or self.max_output_bytes <= 0:
            raise ValueError('Only bounded CPU runners are supported')
        for path, digest in self.file_sha256.items():
            if not Path(path).is_absolute(): raise ValueError('Absolute pinned runner paths required')
            require_sha(digest)
        object.__setattr__(self, 'command', tuple(self.command))
        object.__setattr__(self, 'authorization_refs', tuple(self.authorization_refs))
        object.__setattr__(self, 'file_sha256', MappingProxyType(dict(self.file_sha256)))


class QualifiedRunner:
    def __init__(self, spec, policy): self.spec, self.policy = spec, policy

    @property
    def runner_id(self): return self.spec.runner_id

    def _qualified(self,budget,deadline):
        # Qualification pins executable plus every file argument; shell, module,
        # inline-code and GPU launches cannot sneak through a numeric budget.
        forbidden = {'sh','bash','zsh','fish','env','pip','pip3','conda','ssh','docker','podman'}
        if Path(self.spec.command[0]).name in forbidden: raise IntegrityError('Shell, installer and remote runners are prohibited')
        for arg in self.spec.command:
            if arg.startswith('-'): raise IntegrityError('Runner switches require an explicit wrapper script')
            path = Path(arg)
            if not path.is_absolute() or str(path) not in self.spec.file_sha256:
                raise IntegrityError('Every runner command file must be exactly hash-pinned')
        for path, digest in self.spec.file_sha256.items():
            h=hashlib.sha256()
            with Path(path).open('rb') as stream:
                while True:
                    budget.remaining_seconds()
                    if budget.clock()>=deadline:raise BudgetExceeded('CPU task qualification time exhausted')
                    chunk=stream.read(1024**2)
                    if not chunk:break
                    h.update(chunk)
            if h.hexdigest() != digest:
                raise IntegrityError('Qualified runner file hash mismatch')

    def run(self, request, *, budget):
        def failed(reason, status='failed'):
            return EvidenceDecision(request, True, status, reason, execution='not_run' if status in {'pending_authorization','unavailable'} else 'failed_cpu')
        if (budget.policy != self.policy or self.policy.mode != 'bounded-retrieval'
                or self.spec.runner_id not in self.policy.allowed_runner_ids
                or not set(self.spec.authorization_refs) <= set(self.policy.authorization_refs)):
            return failed('Current applicable runner authorization is required', 'pending_authorization')
        try:
            deadline=budget.clock()+self.policy.task_seconds
            self._qualified(budget,deadline)
            with budget.cpu_runner(self.spec.runner_id) as seconds:
                remaining=min(seconds,deadline-budget.clock())
                if remaining<=0:raise BudgetExceeded('CPU task time exhausted during qualification')
                return self._execute(request, remaining, budget)
        except FileNotFoundError as exc:
            return failed('Qualified runner dependency unavailable: '+str(exc),'unavailable')
        except (IntegrityError, BudgetExceeded, OSError, ValueError) as exc:
            return failed(str(exc))

    def _execute(self, request, seconds, budget):
        deadline = budget.clock() + seconds
        # A small clean environment blocks inherited credential and GPU settings.
        env = {'PATH': '/usr/bin:/bin', 'LANG': 'C.UTF-8', 'CUDA_VISIBLE_DEVICES': '',
               'HIP_VISIBLE_DEVICES': '', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
               'OPENBLAS_NUM_THREADS': '1', 'PYTHONNOUSERSITE': '1'}
        proc = subprocess.Popen(self.spec.command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, env=env, start_new_session=True)
        selector = selectors.DefaultSelector()
        output, errors = bytearray(), bytearray()
        try:
            payload = canonical_bytes(request)
            if len(payload) > 1024**2: raise IntegrityError('Runner request exceeds input byte budget')
            # Request bodies are small bounded JSON, but pipe writes must also
            # respect the task deadline if the qualified process never reads.
            os.set_blocking(proc.stdin.fileno(), False)
            pending = memoryview(payload)
            selector.register(proc.stdin, selectors.EVENT_WRITE, 'stdin')
            for stream, name in ((proc.stdout,'stdout'),(proc.stderr,'stderr')):
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ, name)
            while selector.get_map():
                remaining = min(deadline-budget.clock(), budget.remaining_seconds())
                if remaining <= 0: raise BudgetExceeded('CPU task time budget exhausted')
                for key, _ in selector.select(min(.05, remaining)):
                    stream = key.fileobj
                    if key.data == 'stdin':
                        try: written = os.write(stream.fileno(), pending)
                        except BrokenPipeError: written = len(pending)
                        pending = pending[written:]
                        if not pending:
                            selector.unregister(stream); stream.close()
                    else:
                        chunk = os.read(stream.fileno(), min(65536, self.spec.max_output_bytes-len(output)-len(errors)+1))
                        if not chunk:
                            selector.unregister(stream); stream.close()
                        else:
                            (output if key.data == 'stdout' else errors).extend(chunk)
                            if len(output)+len(errors) > self.spec.max_output_bytes:
                                raise BudgetExceeded('CPU output byte budget exhausted')
            remaining = min(deadline-budget.clock(), budget.remaining_seconds())
            if remaining <= 0: raise BudgetExceeded('CPU task time budget exhausted')
            try: code = proc.wait(timeout=remaining)
            except subprocess.TimeoutExpired as exc: raise BudgetExceeded('CPU task time budget exhausted') from exc
            if code: raise IntegrityError(f'Qualified CPU runner exited with status {code}')
            try: result = decision_from_dict(json.loads(output))
            except (ValueError, TypeError, KeyError) as exc: raise IntegrityError('Invalid typed CPU result') from exc
            if result.request != request or result.qualification != 'diagnostic' or result.scoring_used:
                raise IntegrityError('CPU runner result identity or qualification mismatch')
            return replace(result, execution='cpu:' + self.spec.runner_id,
                           provenance={**result.provenance, 'runner_id':self.spec.runner_id},
                           metadata={**result.metadata, 'authorization_refs':self.spec.authorization_refs,
                                     'runner_file_sha256':dict(self.spec.file_sha256)})
        finally:
            selector.close()
            # Descendants may retain pipes after the immediate process exits.
            # Always terminate the dedicated process group, including that case.
            try: os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            proc.wait()
            for stream in (proc.stdin,proc.stdout,proc.stderr):
                if stream is not None and not stream.closed: stream.close()
