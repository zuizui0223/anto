"""Read only official public OSF metadata/download links, preserving v26 hashes.

A 401/403 is not bypassed. No alternate data, mirror, credential or raw-source
redistribution. The local cache contains only exact size/MD5/SHA256 matches.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import re
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

PROJECT='zg9cu'


def metadata_download(document, record):
    try:
        data=document['data']; attrs=data['attributes']
        node=data['relationships']['node']['links']['related']['href']
        link=data['links']['download']
        if (data['id']!=record['file_id'] or data['type']!='files' or
                attrs['kind']!='file' or attrs['provider']!='osfstorage' or
                attrs['size']!=record['size'] or
                attrs['materialized_path']!=record['source_path']):
            raise ValueError('OSF_METADATA_SOURCE_IDENTITY_MISMATCH')
        n=urlsplit(node); u=urlsplit(link)
        if (n.scheme!='https' or n.hostname!='api.osf.io' or
                n.path!=f'/v2/nodes/{PROJECT}/' or n.username or n.password):
            raise ValueError('OSF_METADATA_PROJECT_MISMATCH')
        allowed=(u.hostname=='files.osf.io' and u.path.rstrip('/')==
                 f"/v1/resources/{PROJECT}/providers/osfstorage/{record['file_id']}") or (
                 u.hostname=='osf.io' and (
                     re.fullmatch(r'/download/[A-Za-z0-9]{5,24}/?',u.path) or
                     re.fullmatch(r'/[A-Za-z0-9]{5}/download/?',u.path)))
        if u.scheme!='https' or u.username or u.password or not allowed:
            raise ValueError('OSF_UNREGISTERED_DOWNLOAD_LINK')
        return link
    except (KeyError,TypeError) as exc:
        raise ValueError('OSF_PUBLIC_METADATA_UNAVAILABLE_OR_WRONG_SCHEMA') from exc


def retrieve(url, maximum):
    req=Request(url,headers={'User-Agent':'ANTO-v31-public-source-retrieval/1.0',
                             'Accept':'application/json,text/plain,text/csv,*/*;q=0.5'})
    with urlopen(req,timeout=60) as response:
        if urlsplit(response.geturl()).scheme!='https':
            raise ValueError('SOURCE_NON_HTTPS_REDIRECT')
        body=response.read(maximum+1)
    if not body or len(body)>maximum:
        raise ValueError('SOURCE_RESPONSE_SIZE')
    return body


def main():
    from audit_petunieae_2023_processed_files_v26 import FILES,verify_original
    p=argparse.ArgumentParser()
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--status',type=Path,required=True)
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    records={k:v for k,v in FILES.items() if k!='pigment_cluster_source'}
    def one(item):
        name,record=item
        status={'file_id':record['file_id'],'expected_bytes':record['size']}
        try:
            raw=retrieve(f"https://api.osf.io/v2/files/{record['file_id']}/",1000000)
            doc=json.loads(raw)
            status['api_download_host']=urlsplit(doc.get('data',{}).get('links',{}).get('download','')).hostname
            status['api_download_path']=urlsplit(doc.get('data',{}).get('links',{}).get('download','')).path
            link=metadata_download(doc,record)
            status['metadata_verified']=True
            status['official_download_host']=urlsplit(link).hostname
            body=retrieve(link,1000000)
            verify_original(body,record)
            (args.cache/(record['file_id']+'.raw')).write_bytes(body)
            status['status']='ORIGINAL_BYTES_SIZE_MD5_SHA256_VERIFIED'
        except Exception as exc:
            status['status']='TECHNICAL_SOURCE_NOT_ADMITTED'
            status['error_type']=type(exc).__name__
            status['http_status']=getattr(exc,'code',None)
            if isinstance(exc,ValueError):status['reason']=str(exc)
        return name,status
    with ThreadPoolExecutor(max_workers=4) as pool:
        statuses=dict(pool.map(one,records.items()))
    result={'source':'official_public_OSF_API_returned_download_links',
            'no_permission_bypass_or_substitute_data':True,'files':statuses,
            'all_originals_verified':all(x['status']=='ORIGINAL_BYTES_SIZE_MD5_SHA256_VERIFIED'
                                         for x in statuses.values())}
    args.status.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('ANTO_V31_OSF_METADATA_TRANSPORT_STATUS',json.dumps(result,sort_keys=True))
    if not result['all_originals_verified']:
        raise SystemExit('Original public source unavailable: no empirical v31 result generated.')

if __name__=='__main__':main()
