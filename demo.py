"""Run the real routing CLI with an explicitly synthetic model bundle."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from demo_fixture import fixture


def run_demo():
    with TemporaryDirectory() as temp:
        root = Path(temp)
        bundle = root / "baseline"
        bundle.mkdir()
        _, pin = fixture(bundle)
        request = root / "request.json"
        request.write_text(json.dumps({
            "chromosome": "1", "position_1based": 10, "ref": "A", "alt": "G"
        }))
        config = root / "config.json"
        config.write_text(json.dumps({
            "schema": "routed-missense-config/v1",
            "baseline": {"bundle": str(bundle), "sha256": pin},
            "cache": str(root / "cache"),
        }))
        result = subprocess.run(
            [sys.executable, "-B", "-m", "evee_next.biology.routed", str(request),
             "--config", str(config), "--policy", "cache-only"],
            check=True, capture_output=True, text=True, timeout=30,
        )
        return json.loads(result.stdout)


if __name__ == "__main__":
    print(json.dumps(run_demo(), indent=2))
