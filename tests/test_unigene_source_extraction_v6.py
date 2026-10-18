"""Offline regression of 2020 Cirsium supplement structure and annotation limits."""
from io import BytesIO
from pathlib import Path
import sys
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from inspect_cirsium_japonicum_unigene_annotation_v6 import (
    fasta_report,worksheet_report,inspect_public_archive,
)

NS="http://schemas.openxmlformats.org/spreadsheetml/2006/main"


def xlsx_fixture():
    shared=[
        "unigene_id","blast_annotation",
        "unigene_1","dihydroflavonol 4-reductase protein",
        "unigene_2","unannotated enzyme"
    ]
    ss=f'<sst xmlns="{NS}" count="{len(shared)}" uniqueCount="{len(shared)}">'+''.join(
        "<si><t>"+v+"</t></si>" for v in shared)+"</sst>"
    table=f'''<worksheet xmlns="{NS}"><sheetData>
        <row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>
        <row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2" t="s"><v>3</v></c></row>
        <row r="3"><c r="A3" t="s"><v>4</v></c><c r="B3" t="s"><v>5</v></c></row>
    </sheetData></worksheet>'''
    xlsx=BytesIO()
    with zipfile.ZipFile(xlsx,"w",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("xl/sharedStrings.xml",ss)
        z.writestr("xl/worksheets/sheet1.xml",table)
        z.writestr("xl/workbook.xml","<workbook/>")
    return xlsx.getvalue()


class SupplementExtractionTest(unittest.TestCase):
    def test_header_count_and_preamble(self):
        b=b"Supporting Information\n>unigene_1\nATGCTGA\n>unigene_2\nCCCTAA\n"
        p=fasta_report(b)
        self.assertEqual(p["n_fasta_headers"],2)
        self.assertEqual(p["total_length_nucleotide_compatible_lines"],13)
        self.assertTrue(p["first_fasta_headers"][0].startswith(">unigene_1"))

    def test_excel_xml_annotations_are_not_orthology(self):
        a=worksheet_report(xlsx_fixture())
        self.assertEqual(a["workbook_row_count"],3)
        self.assertEqual(a["candidate_annotation_row_count_at_least"],1)
        self.assertEqual(a["candidate_annotation_hit_rows_capped"][0]["row"],2)
        self.assertTrue(a["candidate_annotation_hit_rows_capped"][0]["no_validated_orthology"])

    def test_complete_small_archive(self):
        b=BytesIO()
        with zipfile.ZipFile(b,"w",compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr("ao0c04001_si_002.txt",b"Source preamble\n>unigene_1\nATGCTGA\n")
            z.writestr("ao0c04001_si_003.xlsx",xlsx_fixture())
        result=inspect_public_archive(b.getvalue())
        self.assertEqual(result["fasta_evidence"]["n_fasta_headers"],1)
        self.assertEqual(result["workbook_annotation_evidence"]["candidate_annotation_row_count_at_least"],1)
        self.assertFalse(result["genuine_orthologs_identified"])
        self.assertFalse(result["focal_ryukyu_data"])
        self.assertFalse(result["reexpression_confirmed"])

    def test_missing_text_supplement_fails(self):
        b=BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("ao0c04001_si_003.xlsx",xlsx_fixture())
        with self.assertRaises(ValueError):
            inspect_public_archive(b.getvalue())


if __name__=="__main__": unittest.main()
