"""Adversarial tests of original Dryad metadata/size/digest identity."""
import copy
import hashlib
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import fetch_original_Phlox_PdPIRK_Dryad_v39 as src


def fake_meta():
    files=[]
    for i,n in enumerate(src.NAMES):
        raw=b"Pistil_ID,Fruit\nA,1\n" + bytes([65+i])*50
        files.append({
            "path":"exp/"+n,"size":len(raw),"digestType":"md5",
            "digest":hashlib.md5(raw).hexdigest(),
            "_links":{"stash:download":{"href":f"/api/v2/files/{100+i}/download"}},
        })
    return {"_embedded":{"stash:files":files},"total":len(files)}


class OriginalDryadSourceTests(unittest.TestCase):
    def test_only_five_small_registered_files(self):
        self.assertEqual(len(src.NAMES),5)
        self.assertIn("vigs_fruitset.csv",src.NAMES)
        self.assertIn("vigs_pollentubes.csv",src.NAMES)
        self.assertNotIn("DFR_floralscans.zip",src.NAMES)
        self.assertLessEqual(src.MAX_FILE_BYTES,150000)
    def test_version_cannot_be_different_DOI(self):
        d={"identifier":"doi:"+src.DOI,"_links":{"stash:version":{"href":"/api/v2/versions/42"}}}
        self.assertEqual(src.official_version_ref(d),"https://datadryad.org/api/v2/versions/42/files")
        d["identifier"]="doi:10.5061/dryad.WRONG"
        with self.assertRaisesRegex(src.SourceGateError,"DOI_IDENTITY"):
            src.official_version_ref(d)
    def test_external_host_or_large_zip_download_url_prohibited(self):
        for uri in ("https://example.org/api/v2/files/1/download",
                    "http://datadryad.org/api/v2/files/1/download",
                    "https://datadryad.org/api/v2/files/3/download?foo=1",
                    "https://datadryad.org/downloads/files.zip"):
            with self.assertRaises(src.SourceGateError):
                src.admitted_api_url(uri)
    def test_complete_original_source_file_inventory(self):
        a=src.original_file_inventory(fake_meta())
        self.assertEqual(set(a),set(src.NAMES))
        self.assertEqual(a["vigs_fruitset.csv"]["official_file_id"],100)
        self.assertTrue(a["vigs_fruitset.csv"]["download_url"].startswith("https://datadryad.org/api/v2/files/"))
    def test_missing_duplicate_or_paginated_files_stop(self):
        x=fake_meta();x["_embedded"]["stash:files"].pop()
        with self.assertRaisesRegex(src.SourceGateError,"SOURCE_CSV_SET_INCOMPLETE"):
            src.original_file_inventory(x)
        x=fake_meta();x["_embedded"]["stash:files"].append(copy.deepcopy(x["_embedded"]["stash:files"][0]))
        with self.assertRaisesRegex(src.SourceGateError,"DUPLICATE_REGISTERED_FILENAME"):
            src.original_file_inventory(x)
        x=fake_meta();x["total"]=10
        with self.assertRaisesRegex(src.SourceGateError,"FILE_PAGE_INCOMPLETE"):
            src.original_file_inventory(x)
    def test_bad_file_uri_or_size_not_admitted(self):
        for field,value in (("size",10**9),("size",0)):
            x=fake_meta();x["_embedded"]["stash:files"][0][field]=value
            with self.assertRaisesRegex(src.SourceGateError,"DECLARED_SIZE"):
                src.original_file_inventory(x)
        x=fake_meta();x["_embedded"]["stash:files"][0]["_links"]["stash:download"]["href"]="https://someother.com/evil"
        with self.assertRaisesRegex(src.SourceGateError,"ORIGINAL_ID"):
            src.original_file_inventory(x)
    def test_digest_mismatch_and_same_length_tamper(self):
        raw=b"Pistil_ID,Fruit\nA,1\n"+b"A"*50
        record=src.original_file_inventory(fake_meta())["vigs_fruitset.csv"]
        proof=src.admit_original_file(raw,record)
        self.assertTrue(proof["original_provider_digest_matched"])
        changed=raw.replace(b"A",b"B")
        self.assertEqual(len(changed),len(raw))
        with self.assertRaisesRegex(src.SourceGateError,"DIGEST_MISMATCH"):
            src.admit_original_file(changed,record)
    def test_missing_original_digest_not_sha_derived_fake_authentication(self):
        raw=b"Pistil_ID,Fruit\nA,1\n"+b"A"*50
        a=src.original_file_inventory(fake_meta())["vigs_fruitset.csv"]
        a["original_manifest_digest"]=None
        with self.assertRaisesRegex(src.SourceGateError,"NO_ORIGINAL_PROVIDER_DIGEST"):
            src.admit_original_file(raw,a)
    def test_sha256_original_supported(self):
        raw=b"Pistil_ID,Fruit\nA,1\n"+b"A"*50
        a=src.original_file_inventory(fake_meta())["vigs_fruitset.csv"]
        a["original_manifest_digest_type"]="sha256"
        a["original_manifest_digest"]=hashlib.sha256(raw).hexdigest()
        self.assertTrue(src.admit_original_file(raw,a)["original_provider_digest_matched"])


if __name__=="__main__":unittest.main()
