"""Public OSF metadata transport, never a source/permission bypass."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import fetch_verified_osf_v31 as m
except ModuleNotFoundError:
    m=None

class SourceTransport(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m,'OSF adapter not yet implemented')
        self.record={'file_id':'a'*24,'size':3,'source_path':'/f.csv'}
        self.doc={'data':{'id':'a'*24,'type':'files',
            'attributes':{'kind':'file','provider':'osfstorage','size':3,'materialized_path':'/f.csv'},
            'relationships':{'node':{'links':{'related':{'href':'https://api.osf.io/v2/nodes/zg9cu/'}}}},
            'links':{'download':'https://files.osf.io/v1/resources/zg9cu/providers/osfstorage/'+'a'*24}}}

    def test_original_official_download_link(self):
        self.assertEqual(m.metadata_download(self.doc,self.record),self.doc['data']['links']['download'])

    def test_official_metadata_short_guid_download_still_uses_pinned_body(self):
        for link in ('https://osf.io/ab12c/download', 'https://osf.io/download/ab12c/'):
            self.doc['data']['links']['download']=link
            self.assertEqual(m.metadata_download(self.doc,self.record),link)

    def test_foreign_project_or_file_cannot_be_substituted(self):
        self.doc['data']['id']='b'*24
        with self.assertRaises(ValueError):m.metadata_download(self.doc,self.record)
        self.doc['data']['id']='a'*24
        self.doc['data']['relationships']['node']['links']['related']['href']='https://api.osf.io/v2/nodes/other/'
        with self.assertRaises(ValueError):m.metadata_download(self.doc,self.record)

    def test_unregistered_size_or_path_fails_closed(self):
        self.doc['data']['attributes']['size']=4
        with self.assertRaises(ValueError):m.metadata_download(self.doc,self.record)
        self.doc['data']['attributes']['size']=3
        self.doc['data']['attributes']['materialized_path']='/replacement.csv'
        with self.assertRaises(ValueError):m.metadata_download(self.doc,self.record)

    def test_permission_or_foreign_url_is_not_overridden(self):
        with self.assertRaises(ValueError):m.metadata_download({'errors':[{'code':'403'}]},self.record)
        for bad in ['http://files.osf.io/a','https://evil.test/a','https://files.osf.io.evil.test/a',
                    'https://files.osf.io/v1/resources/other/providers/osfstorage/'+'a'*24]:
            self.doc['data']['links']['download']=bad
            with self.assertRaises(ValueError):m.metadata_download(self.doc,self.record)

if __name__=='__main__':unittest.main()
