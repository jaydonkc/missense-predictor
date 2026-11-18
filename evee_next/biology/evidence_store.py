"""Immutable inference artifacts, separate from study-target caches.

An entry is visible only after its payload and manifest have been fsynced and
the complete staging directory atomically renamed. All reads use checked byte
snapshots. Failed attempts are receipts, never negative biological evidence.
"""
from contextlib import contextmanager
from dataclasses import replace
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
import math
import uuid

from .routing import EvidenceDecision, EvidenceRequest, canonical_bytes, public_value, require_sha


class CacheIntegrityError(ValueError):
    """An existing entry is incomplete, inconsistent or no longer immutable."""


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def decision_from_dict(value):
    value = dict(value)
    value['request'] = EvidenceRequest(**value['request'])
    return EvidenceDecision(**value)


class EvidenceStore:
    def __init__(self, root, *, max_artifact_bytes=256 * 1024**2,lock_timeout_seconds=30.):
        self.root = Path(root)
        if type(max_artifact_bytes) is not int or max_artifact_bytes <= 0:
            raise ValueError('Positive artifact read limit required')
        self.max_artifact_bytes = max_artifact_bytes
        if isinstance(lock_timeout_seconds,bool) or not isinstance(lock_timeout_seconds,(int,float)) or not math.isfinite(lock_timeout_seconds) or lock_timeout_seconds<=0:
            raise ValueError('Finite positive cache lock timeout required')
        self.lock_timeout_seconds=lock_timeout_seconds
        for name in ('entries', 'blobs', 'locks', 'attempts', 'staging'):
            path = self.root / name
            path.mkdir(parents=True, exist_ok=True)
            if path.is_symlink(): raise CacheIntegrityError('Cache directories cannot be symlinks')

    @contextmanager
    def _lock(self, key):
        fd = os.open(self.root / 'locks' / key, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
        try:
            deadline=time.monotonic()+self.lock_timeout_seconds
            while True:
                try:
                    fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    remaining=deadline-time.monotonic()
                    if remaining<=0:raise TimeoutError('Cache key lock exceeded bounded wait')
                    time.sleep(min(.01,remaining))
            yield
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)

    def _snapshot(self, path, limit):
        try:
            fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
            with os.fdopen(fd, 'rb') as handle:
                data = handle.read(limit + 1)
            if len(data) > limit: raise CacheIntegrityError('Cache payload exceeds configured snapshot limit')
            return data
        except OSError as exc:
            raise CacheIntegrityError('Missing or unsafe cache payload') from exc

    def _entry(self, request):
        path = self.root / 'entries' / request.key
        if not path.exists(): return None
        if path.is_symlink() or not path.is_dir(): raise CacheIntegrityError('Unsafe cache entry')
        try:
            manifest = json.loads(self._snapshot(path / 'manifest.json', 4 * 1024**2))
            if manifest['schema'] != 'inference-evidence-cache/v1' or manifest['key'] != request.key:
                raise CacheIntegrityError('Cache schema or identity mismatch')
            if _sha(canonical_bytes(manifest['decision'])) != manifest['decision_sha256']:
                raise CacheIntegrityError('Decision hash mismatch')
            result = decision_from_dict(manifest['decision'])
            if result.request != request: raise CacheIntegrityError('Cache request mismatch')
            data = self._snapshot(path / 'artifact.bin', self.max_artifact_bytes)
            if len(data) != manifest['artifact_bytes'] or _sha(data) != manifest['artifact_sha256']:
                raise CacheIntegrityError('Artifact hash or length mismatch')
            return result, data, manifest
        except (KeyError, TypeError, ValueError, UnicodeError) as exc:
            if isinstance(exc, CacheIntegrityError): raise
            raise CacheIntegrityError('Invalid cache manifest') from exc

    def lookup(self, request):
        entry = self._entry(request)
        return None if entry is None else replace(entry[0], cache_disposition='hit')

    def read_artifact(self, request):
        """Return the hash-checked immutable bytes bound to this exact request."""
        entry = self._entry(request)
        return None if entry is None else entry[1]

    def get_artifact(self, sha256):
        """Read a shared raw blob by its pinned content hash, without parsing."""
        require_sha(sha256)
        path = self.root / 'blobs' / sha256
        if not path.exists(): return None
        data = self._snapshot(path, self.max_artifact_bytes)
        if _sha(data) != sha256: raise CacheIntegrityError('Shared artifact hash mismatch')
        return data

    @staticmethod
    def _write(path, data):
        with open(path, 'xb') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())

    @staticmethod
    def _sync_directory(path):
        fd = os.open(path, os.O_RDONLY)
        try: os.fsync(fd)
        finally: os.close(fd)

    def _publish(self, request, artifact_bytes, decision):
        if not isinstance(artifact_bytes, bytes): raise TypeError('Immutable artifact bytes required')
        if decision.request != request: raise ValueError('Decision differs from requested identity')
        if len(artifact_bytes) > self.max_artifact_bytes: raise CacheIntegrityError('Artifact too large')
        # Missing/unavailable/failed requests may be retried after prerequisites change.
        # Preserve their receipts without making them valid cache hits.
        if decision.status not in {'available', 'partial', 'not_applicable'}:
            target = self.root / 'attempts' / (request.key + '-' + uuid.uuid4().hex + '.json')
            self._write(target, canonical_bytes(decision))
            return replace(decision, cache_disposition='attempt_recorded')
        stored = replace(decision, cache_disposition='published')
        existing = self._entry(request)
        if existing is not None:
            normalized = replace(existing[0], cache_disposition='published')
            if normalized != stored or existing[1] != artifact_bytes:
                raise CacheIntegrityError('Immutable cache publication conflict')
            return replace(existing[0], cache_disposition='hit')
        digest = _sha(artifact_bytes)
        stage = Path(tempfile.mkdtemp(prefix=request.key + '-', dir=self.root / 'staging'))
        try:
            self._write(stage / 'artifact.bin', artifact_bytes)
            manifest = dict(schema='inference-evidence-cache/v1', key=request.key,
                            artifact_sha256=digest, artifact_bytes=len(artifact_bytes),
                            decision=public_value(stored), decision_sha256=_sha(canonical_bytes(stored)))
            self._write(stage / 'manifest.json', canonical_bytes(manifest))
            self._sync_directory(stage)
            # A blob link is created only from a fully written file; competing
            # requests can safely publish identical content simultaneously.
            blob = self.root / 'blobs' / digest
            try: os.link(stage / 'artifact.bin', blob)
            except FileExistsError:
                if self.get_artifact(digest) != artifact_bytes: raise CacheIntegrityError('Shared blob conflict')
            self._sync_directory(self.root / 'blobs')
            os.rename(stage, self.root / 'entries' / request.key)
            self._sync_directory(self.root / 'entries')
        finally:
            if stage.exists(): shutil.rmtree(stage)
        return stored

    def publish(self, request, artifact_bytes, decision):
        with self._lock(request.key):
            return self._publish(request, artifact_bytes, decision)

    def get_or_compute(self, request, producer):
        """Lock/recheck/compute/publish once for concurrent identical requests."""
        with self._lock(request.key):
            cached = self.lookup(request)
            if cached is not None: return cached
            artifact_bytes, decision = producer()
            return self._publish(request, artifact_bytes, decision)
