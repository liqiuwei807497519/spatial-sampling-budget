"""Evaluate the five named suspension matrices once each, including ineligible barcodes in strict draws."""
from pathlib import Path
from zipfile import ZipFile
import gzip,io,argparse,math,hashlib
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.special import gammaln

def no_event_mean(q,m):
    # Elementary symmetric means with the sequential finite-population update.
    e=np.zeros(m+1,dtype=np.longdouble);e[0]=1
    for n,x in enumerate(q,1):
        k=np.arange(1,min(m,n)+1)
        e[k]=e[k]*(1-k/n)+e[k-1]*x*(k/n)
    return float(e[m])
def score_inputs(inputs,ref,out):
    from contextlib import ExitStack
    rows=[];counts=[]
    if len({s for s,_ in inputs})!=len(inputs):raise ValueError('Duplicate input sample')
    with ExitStack() as stack:
        archives={p:stack.enter_context(ZipFile(p)) for p in set(p for _,p in inputs)}
        features={};common=set(ref.index)
        for sample,path in inputs:
            z=archives[path];entry=f'{sample}/matrix/matrix.mtx.gz'
            folder=entry.rsplit('/',1)[0]+'/'
            ft=pd.read_csv(io.BytesIO(gzip.decompress(z.read(folder+'features.tsv.gz'))),sep='\t',header=None,dtype=str,keep_default_na=False)
            ft.columns=['id','gene','type'][:ft.shape[1]]
            ft.id=ft.id.str.replace(r'\.\d+$','',regex=True)
            if 'type' not in ft:ft['type']='Gene Expression'
            ok=(ft.type=='Gene Expression')&~ft.id.duplicated(keep=False)&(ft.gene==ft.id.map(ref))
            common&=set(ft.loc[ok,'id']);features[entry]=ft
        assert {'ENSG00000157766','ENSG00000139219'}<=common or {'ACAN','COL2A1'}<=set(ref.loc[list(common)])
        for sample,path in inputs:
            z=archives[path];entry=f'{sample}/matrix/matrix.mtx.gz';ft=features[entry]
            sample=entry.split('/')[0]
            folder=entry.rsplit('/',1)[0]+'/'
            bc=gzip.decompress(z.read(folder+'barcodes.tsv.gz')).splitlines()
            assert len(set(bc))==len(bc)
            valid=ft.id.isin(common).to_numpy()
            A=(ft.gene=='ACAN').to_numpy()&valid;B=(ft.gene=='COL2A1').to_numpy()&valid
            assert A.sum()==B.sum()==1
            with z.open(entry) as f,gzip.GzipFile(fileobj=f) as g:
                header=g.readline().lower().strip();assert b'matrix coordinate integer' in header
                line=g.readline()
                while line.startswith(b'%'):line=g.readline()
                nf,nb,nnz=map(int,line.split());assert nf==len(ft) and nb==len(bc)
                sums=np.zeros((3,nb),dtype=np.int64);seen=0
                for chunk in pd.read_csv(g,sep=r'\s+',header=None,chunksize=400000,dtype=np.int64):
                    ri,ci,val=chunk.to_numpy().T;ri=ri-1;ci=ci-1
                    assert ri.min()>=0 and ri.max()<nf and ci.min()>=0 and ci.max()<nb and val.min()>0
                    for j,mask in enumerate((valid,A,B)):
                        use=mask[ri];sums[j]+=np.bincount(ci[use],weights=val[use],minlength=nb).astype(np.int64)
                    seen+=len(val)
                assert seen==nnz
            L,a,b=sums;k=521;m=12;ok=(L-a-b>=k)
            passp=math.prod((ok.sum()-j)/(nb-j) for j in range(m)) if ok.sum()>=m else 0.
            if passp:
                miss=[hypergeom.pmf(0,L[ok],g[ok],k) for g in [a,b,a+b]]
                both=1-no_event_mean(miss[0],m)-no_event_mean(miss[1],m)+no_event_mean(miss[2],m)
                single=np.clip(1-miss[0]-miss[1]+miss[2],0,1)
                within=1-no_event_mean(1-single,m)
            else:both=within=np.nan
            rows.append(dict(library=sample,common_gene_axis=len(common),released=nb,eligible=int(ok.sum()),twelve_pass=passp,dual_among_full=both,strict_dual=passp*both if passp else 0.,within_barcode_strict=passp*within if passp else 0.,metadata_region='disc tissue, compartment unspecified'))
            counts.append(pd.DataFrame(dict(library=sample,total=L,ACAN=a,COL2A1=b)))
            print(rows[-1],flush=True)
    out.mkdir(exist_ok=True,parents=True)
    pd.DataFrame(rows).to_csv(out/'OMIX017533_suspension.csv',index=False)
    pd.concat(counts).to_csv(out/'OMIX017533_barcode_counts.csv.gz',index=False)

def main():
    from ngdc_files import mapped_inputs
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path,help='Directory containing the five official archives, or a checksum-verified legacy combined archive')
    ap.add_argument('--download',action='store_true',help='Download any missing official file and check every SHA-256')
    ap.add_argument('--axis',type=Path,default=Path(__file__).parent/'reference_21479_to_20512.csv')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    ref=pd.read_csv(args.axis).query('eligible==1').set_index('stable_gene_id').symbol
    if len(ref)!=20512 or ref.index.has_duplicates:raise ValueError('Invalid frozen reference')
    score_inputs(mapped_inputs(args.source,args.download),ref,args.out)
if __name__=='__main__':main()
