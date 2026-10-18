"""Replay published 2020 Cirsium transcript protein-similarity triage offline.

No new RNA-seq downloads. This is NOT Cirsium brevicaule/irumtiense genotyping,
ortholog certification, anthocyanin expression or phenotypic re-evolution.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_2020_homology_triage_v6 import run_study


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--transcripts",type=Path,
                   default=ROOT/"data/molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta")
    p.add_argument("--metadata",type=Path,
                   default=ROOT/"data/molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json")
    p.add_argument("--reference",type=Path,
                   default=ROOT/"data/molecular/asteraceae_reference_cds_v5.fasta")
    p.add_argument("--output",type=Path,default=None)
    args=p.parse_args()
    result=run_study(args.transcripts,args.metadata,args.reference)
    content=json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")
    else:
        print(content)
    return result


if __name__=="__main__":
    main()
