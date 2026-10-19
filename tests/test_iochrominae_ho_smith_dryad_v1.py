"""Synthetic negative controls for an official Dryad source-only intake."""
import io
import json
import hashlib
import urllib.error
from unittest.mock import patch
from pathlib import Path
import sys
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_iochrominae_ho_smith_original_dryad_v1 as d

def fixture(dfr_header="28 1179",missing=None):
    f={"EVOL_CHI_SEQFILE.txt":"28 717\n"+"ATGC"*350,
       "EVOL_DFR_SEQFILE.txt":dfr_header+"\n"+"ATGC"*350,
       "F3H_SEQFILE.txt":"28 1113\n"+"ATGC"*350,
       "README_for_EVOL_CHI_SEQFILE.txt":"original source legend"}
    if missing:f.pop(missing)
    bio=io.BytesIO()
    with zipfile.ZipFile(bio,"w") as z:
        for p,s in f.items():z.writestr(p,s)
    return bio.getvalue()


class SourceIntakeTests(unittest.TestCase):
    def test_exact_three_original_gene_file_names_and_headers(self):
        report=d.inspect_zip(fixture())
        self.assertEqual(set(report["target_gene_file_metadata"]),set(d.EXPECTED_FILES))
        self.assertEqual(report["target_gene_file_metadata"]["EVOL_CHI_SEQFILE.txt"]["declared_alignment_sites"],717)
        self.assertEqual(report["target_gene_file_metadata"]["EVOL_DFR_SEQFILE.txt"]["declared_sequences"],28)
        self.assertEqual(report["target_gene_file_metadata"]["F3H_SEQFILE.txt"]["declared_alignment_sites"],1113)

    def test_missing_one_gene_alignment_is_not_success(self):
        with self.assertRaisesRegex(ValueError,"missing"):d.inspect_zip(fixture(missing="F3H_SEQFILE.txt"))

    def test_source_identity_and_limits(self):
        m={"identifier":"doi:10.5061/dryad.49rs8","versionStatus":"published"}
        r=d.evaluate(m,fixture())
        self.assertEqual(r["biological_inference_level"],"SOURCE_FILE_RETRIEVAL_AND_HEADER_ONLY")
        self.assertEqual(r["outcomes_causally_measured_here"],[])
        self.assertFalse(r["original_sequences_rehosted"])
        self.assertTrue(r["no_ancestral_regain_or_phenotypic_persistence_inferred"])

    def test_wrong_doi_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"dataset"):d.evaluate({"identifier":"doi:10.5061/dryad.WRONG"},fixture())

    def test_url_restriction(self):
        with self.assertRaisesRegex(ValueError,"official"):d.get_bytes("https://example.com/other.txt")

    def test_public_file_api_fallback_reconstructs_gene_files_with_source_digests(self):
        from hashlib import sha256
        with zipfile.ZipFile(io.BytesIO(fixture())) as z:
            source={k:z.read(k) for k in d.EXPECTED_FILES}
        metadata={"identifier":"doi:10.5061/dryad.49rs8",
                  "versionStatus":"published",
                  "_links":{"stash:version":{"href":"/api/v2/versions/12345"}}}
        records=[{"path":k,"digest":sha256(v).hexdigest(),
                  "digestType":"sha-256",
                  "_links":{"stash:download":{"href":f"/api/v2/files/{i+1000}/download"}}}
                 for i,(k,v) in enumerate(source.items())]
        def mocked(url):
            if url==d.DOWNLOAD_URL:
                raise urllib.error.HTTPError(url,401,"bundle auth required",{},None)
            if url=="https://datadryad.org/api/v2/versions/12345/files":
                return json.dumps({"total":3,"_embedded":{"stash:files":records}}).encode()
            for i,(_,v) in enumerate(source.items()):
                if url==f"https://datadryad.org/api/v2/files/{i+1000}/download":
                    return v
            raise AssertionError(url)
        with patch.object(d,"get_bytes",side_effect=mocked):
            b,mode=d.public_original_material(metadata)
        self.assertEqual(mode,"official_registered_file_API_fallback")
        r=d.evaluate(metadata,b,mode)
        self.assertIsNone(r["original_zip_sha256"])
        self.assertEqual(set(r["source_files"]["target_gene_file_metadata"]),set(d.EXPECTED_FILES))

    def test_file_api_digest_mismatch_abstains(self):
        metadata={"_links":{"stash:version":{"href":"/api/v2/versions/12345"}}}
        row={"path":"EVOL_CHI_SEQFILE.txt","digest":"0"*64,"digestType":"sha-256",
             "_links":{"stash:download":{"href":"/api/v2/files/1000/download"}}}
        def mocked(url):
            if url.endswith("/files"): return json.dumps({"total":1,"_embedded":{"stash:files":[row]}}).encode()
            return b"A"*2000
        with patch.object(d,"get_bytes",side_effect=mocked):
            with self.assertRaisesRegex(ValueError,"SHA256 mismatch"):
                d.download_registered_files(metadata)

    def test_access_denial_records_explicit_metadata_only_abstention(self):
        meta={"identifier":"doi:10.5061/dryad.49rs8","title":"Original Ho-Smith",
              "_links":{"stash:version":{"href":"/api/v2/versions/12345"}}}
        names=list(d.EXPECTED_FILES)
        files=[{"path":name,"size":12345,"digest":"abc"+str(i),"digestType":"sha-256",
                "_links":{"stash:download":{"href":f"/api/v2/files/{i+2222}/download"}}}
               for i,name in enumerate(names)]
        def mocked(url):
            if url=="https://datadryad.org/api/v2/versions/12345/files":
                return json.dumps({"total":3,"_embedded":{"stash:files":files}}).encode()
            raise AssertionError(url)
        with patch.object(d,"get_bytes",side_effect=mocked):
            receipt=d.blocked_source_receipt(meta,"SOURCE_FILE_ACCESS_RESTRICTED_HTTP_401")
        self.assertEqual(receipt["scientific_admission"],"SOURCE_ONLY_METADATA_FILES_UNAVAILABLE")
        self.assertFalse(receipt["source_bytes_downloaded"])
        self.assertFalse(receipt["per_file_sha256_computed"])
        self.assertFalse(receipt["coding_integrity_independently_replayed"])
        self.assertFalse(receipt["phenotypic_regain_or_fitness_estimated"])
        self.assertEqual(set(receipt["original_files_metadata_only"]),set(d.EXPECTED_FILES))

    def test_missing_version_link_is_not_caught_with_fake_file(self):
        with self.assertRaisesRegex(RuntimeError,"VERSION_METADATA"):
            d.download_registered_files({})

    def test_unknown_format_must_remain_unknown(self):
        x=d.read_alignment_header(b"not-an-alignment\nACGT\n")
        self.assertEqual(x["format"],"UNCLASSIFIED")

    def test_fasta_is_not_silently_interpreted_as_coding(self):
        x=d.read_alignment_header(b">species_A\nATGC\n>species_B\nATGC\n")
        self.assertEqual(x["format"],"FASTA")
        self.assertIsNone(x["declared_alignment_sites"])

if __name__=="__main__":unittest.main()
