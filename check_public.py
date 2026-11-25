"""Verify both strict endpoints independently using Newton identities on the barcode totals."""
import argparse,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

def elementary(q,m):
    if len(q)<m:raise ValueError("Too few eligible barcodes")
    x=np.asarray(q,dtype=np.longdouble);n=len(x)
    powers=[0]+[np.sum(x**k)/np.longdouble(n)**k for k in range(1,m+1)]
    e=[np.longdouble(1.)]
    for k in range(1,m+1):
        e.append(sum((-1)**(j-1)*e[k-j]*powers[j] for j in range(1,k+1))/k)
    return float(e[m]/(np.longdouble(math.comb(n,m))/np.longdouble(n)**m))
def verify(folder,expected=None):
    raw=pd.read_csv(folder/'OMIX017533_barcode_counts.csv.gz')
    scores=pd.read_csv(folder/'OMIX017533_suspension.csv').set_index('library')
    assert set(scores.index)==set(raw.library)=={f'IVD{i}' for i in range(1,6)}
    for name,z in raw.groupby('library'):
        f=z[z.total-z.ACAN-z.COL2A1>=521]
        passp=math.comb(len(f),12)/math.comb(len(z),12)
        q=[hypergeom.pmf(0,f.total,g,521) for g in (f.ACAN,f.COL2A1,f.ACAN+f.COL2A1)]
        both=1-elementary(q[0],12)-elementary(q[1],12)+elementary(q[2],12)
        single=np.clip(1-q[0]-q[1]+q[2],0,1)
        within=passp*(1-elementary(1-single,12))
        np.testing.assert_allclose([passp,both,passp*both,within],scores.loc[name,['twelve_pass','dual_among_full','strict_dual','within_barcode_strict']].to_numpy(float),atol=1e-12,rtol=0)
        assert within<=passp*both+1e-12
    if expected:
        ref=pd.read_csv(expected).set_index('library').sort_index()
        pd.testing.assert_frame_equal(scores.sort_index(),ref,atol=1e-12,rtol=1e-12)
    print('Verified both endpoints in five libraries and',len(raw),'released barcodes.')
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output',type=Path);ap.add_argument('--expected',type=Path,default=None)
    args=ap.parse_args();verify(args.output,args.expected)
