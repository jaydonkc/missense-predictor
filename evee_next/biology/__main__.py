import argparse
import json
from pathlib import Path
from .runtime import CachedEffectProvider


def main():
    p=argparse.ArgumentParser(description='Inspect explicit chemistry/context/abundance stages; no clinical prediction from a protein-only request')
    p.add_argument('gene',choices=['ASPA','G6PD']);p.add_argument('position',type=int)
    p.add_argument('ref');p.add_argument('alt');p.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[2])
    a=p.parse_args()
    report=CachedEffectProvider(a.project_root).query(a.gene,a.position,a.ref,a.alt)
    print(json.dumps(report,sort_keys=True,indent=2,allow_nan=False))


if __name__=='__main__':main()
