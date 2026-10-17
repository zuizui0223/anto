"""Offline replay of independently retrieved PRJNA1311153 public *metadata*.

python scripts/run_public_sra_intake_v4.py --out results/public_sra_audit_v4.json

No raw reads are opened; no external API is called in routine CI. A future
source refresh should be staged and compared by immutable run/sample keys,
rather than silently replacing a frozen source snapshot.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from anto_sra_intake import audit


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=None)
    p.add_argument("--manifest",type=Path,default=ROOT/"data/public_sra/prjna1311153_run_manifest_20261017.csv")
    p.add_argument("--morph",type=Path,default=ROOT/"data/public_sra/takaoense_aza3_morph_crosswalk_20261017.csv")
    p.add_argument("--platform",type=Path,default=ROOT/"data/public_sra/prjna1311153_ena_platform_fingerprints_20261017.json")
    args=p.parse_args()
    result=audit(args.manifest,args.morph,args.platform)
    result["analysis_version"]="v4"
    result["analysis_type"]="public_metadata_only"
    result["registry_snapshot_date"]="2026-10-17"
    s=json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(s,encoding="utf-8")
    else:print(s)
    return result


if __name__=="__main__":
    main()
