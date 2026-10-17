"""Offline v5: check published Asteraceae reference CDS and focal assay readiness.

No downloads or analysis of the actual six focal RNA-seq samples occurs here.
Optional --candidates accepts externally generated candidate CDSs and performs
sequence-level QC only, never genuine allele-function or ancestral-state inference.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_molecular_screen_v5 import screen


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--reference-fasta",type=Path,default=ROOT/"data/molecular/asteraceae_reference_cds_v5.fasta")
    parser.add_argument("--reference-provenance",type=Path,default=ROOT/"data/molecular/asteraceae_reference_provenance_v5.json")
    parser.add_argument("--register",type=Path,default=ROOT/"data/molecular/ryukyu_gene_screen_register_v5.csv")
    parser.add_argument("--run-manifest",type=Path,default=ROOT/"data/public_sra/prjna1311153_run_manifest_20261017.csv")
    parser.add_argument("--candidates",type=Path,default=None,
                        help="Optional proposed CDS FASTA, headers >SRR|locus|contig; no biological orthology implied")
    parser.add_argument("--output",type=Path,default=None)
    args=parser.parse_args()
    record=screen(args.reference_fasta,args.reference_provenance,
                  args.register,args.run_manifest,args.candidates)
    text=json.dumps(record,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(exist_ok=True,parents=True)
        args.output.write_text(text,encoding="utf-8")
    else:print(text)
    return record


if __name__=="__main__":main()
