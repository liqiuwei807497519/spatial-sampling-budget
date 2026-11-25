"""Pinned registry and resumable HTTPS downloads for the five official OMIX files."""
import csv
import hashlib
import os
import shutil
import tempfile
import subprocess
from pathlib import Path
from urllib.request import Request,urlopen
from concurrent.futures import ThreadPoolExecutor
from zipfile import ZipFile

MANIFEST=Path(__file__).parent/'OMIX017533_files.csv'
def records():
    with MANIFEST.open(encoding='utf8') as stream:return list(csv.DictReader(stream))
def sha256(p):
    with Path(p).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def verify(p,row):
    if not p.is_file() or p.stat().st_size!=int(row['bytes']) or sha256(p)!=row['sha256']:
        raise ValueError(f"File integrity check failed for {p.name}; use the matching official archive.")
def acquire(row,folder):
    p=folder/row['file']
    if p.exists():verify(p,row);return p
    size=int(row['bytes']);chunk=4*1024**2
    parts=folder/('.'+p.name+'.parts');parts.mkdir(exist_ok=True)
    def read_part(a):
        b=min(size,a+chunk)-1;part=parts/str(a)
        if part.is_file() and part.stat().st_size==b-a+1:return part
        for attempt in range(4):
            try:
                if shutil.which('curl') and attempt in (0,2,3):
                    with tempfile.TemporaryDirectory(dir=parts) as wd:
                        h=Path(wd)/'headers';payload=Path(wd)/'payload'
                        r=subprocess.run(['curl','-fsSL','--connect-timeout','12','--max-time','40','-r',f'{a}-{b}','-D',str(h),'-o',str(payload),row['url']+f'?jtm={a}&v={attempt}'],capture_output=True)
                        if r.returncode or f'bytes {a}-{b}/{size}' not in h.read_text(errors='ignore').lower() or payload.stat().st_size!=b-a+1:
                            raise IOError('Incomplete HTTPS range')
                        os.replace(payload,part)
                    return part
                req=Request(row['url']+f'?part={a}&try={attempt}',headers={'Range':f'bytes={a}-{b}'})
                with urlopen(req,timeout=50) as resp:
                    if resp.getheader('Content-Range','').lower()!=f'bytes {a}-{b}/{size}':raise IOError('Unexpected HTTP byte range')
                    data=resp.read(b-a+2)
                if len(data)!=b-a+1:raise IOError('Short download')
                with tempfile.NamedTemporaryFile(dir=parts,delete=False) as tmp:tmp.write(data)
                os.replace(tmp.name,part)
                return part
            except (OSError,ValueError):
                if attempt==3:raise
    with ThreadPoolExecutor(6) as pool:segments=list(pool.map(read_part,range(0,size,chunk)))
    with tempfile.NamedTemporaryFile(dir=folder,delete=False) as tmp:
        for part in segments:
            with part.open('rb') as reader:shutil.copyfileobj(reader,tmp)
    temp=Path(tmp.name)
    try:
        verify(temp,row);os.replace(temp,p)
    finally:temp.unlink(missing_ok=True)
    return p
def mapped_inputs(source,download=False):
    source=Path(source); table=records()
    if source.is_file():
        row=next(r for r in table if r['library']=='IVD1')
        verify(source,row)
        return [(r['library'],source) for r in table]
    if download:source.mkdir(exist_ok=True,parents=True)
    mapped=[]
    for row in table:
        p=acquire(row,source) if download else source/row['file']
        verify(p,row)
        mapped.append((row['library'],p))
    if len({sample for sample,_ in mapped})!=len(table):raise ValueError('Duplicate sample assignments')
    return sorted(mapped)
