"""Compare finite-frame designs at a fixed total RNA budget. Run locally on the supplied rat data."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.spatial.distance import cdist
from sampling_design import POLICIES,draws,partition,benchmark

GENES=['Acan','Col2a1','Col1a1','Col9a1','Comp','Dcn','Bgn','Fmod','Sox9','Krt19','Mmp3','Mmp13','Adamts5']
PLANS=[(6,100),(12,50),(20,30),(30,20)]
METHODS=(*POLICIES,'maximin','inhibition')
BASE=Path(__file__).resolve().parent

def farthest(frame,m,n,seed):
    """Random first point followed by greedy farthest point with stable index tie-breaking."""
    xy=frame[['x','y']].to_numpy(float)
    rng=np.random.default_rng(seed)
    starts=rng.integers(len(xy),size=n)
    unique,back=np.unique(starts,return_inverse=True)
    out=np.empty((len(unique),m),dtype=int)
    dd=cdist(xy,xy,'sqeuclidean')
    for row,first in zip(out,unique):
        row[0]=first
        nearest=dd[first].copy()
        for j in range(1,m):
            chosen=int(np.argmax(nearest))
            row[j]=chosen
            np.minimum(nearest,dd[chosen],out=nearest)
    sel=out[back]
    dark=frame.stratum.to_numpy()[sel]=='dark'
    a=dark.sum(1); b=m-a
    w=dark/np.maximum(a[:,None],1)-(~dark)/np.maximum(b[:,None],1)
    return sel,w,(a*b>0)

def policies(f,m,n,seed,method):
    if method=='maximin':return farthest(f,m,n,seed)
    if method=='inhibition':
        # Fixed-radius sequential inhibition on a finite set with randomized visitation order.
        xy=f[['x','y']].to_numpy()
        ds=cdist(xy,xy)
        np.fill_diagonal(ds,np.inf)
        pitch=np.median(ds.min(axis=1))
        radius=.55*pitch*np.sqrt(len(xy)/m)
        neighbours=[np.flatnonzero(row<radius) for row in ds]
        rng=np.random.default_rng(seed)
        indices=np.zeros((n,m),int)
        for row in indices:
            excluded=np.zeros(len(xy),bool)
            filled=0
            for j in rng.permutation(len(xy)):
                if excluded[j]:continue
                row[filled]=j;filled+=1
                if filled==m:break
                excluded[neighbours[j]]=True
            if filled!=m:raise ValueError('Inhibition radius cannot support requested layout')
        flag=f.stratum.to_numpy()[indices]=='dark'
        nd=flag.sum(1);nl=m-nd
        w=flag/np.maximum(nd,1)[:,None]-(~flag)/np.maximum(nl,1)[:,None]
        return indices,w,nd*nl>0
    return draws(f,m,n,seed,method)

def evaluate(f,sel,w,ok,k,pairs):
    total=f.total.to_numpy()
    signs=f.stratum.eq('dark').to_numpy()
    result=[]
    for group,(g,h) in pairs.items():
        a,b=f[g].to_numpy(),f[h].to_numpy()
        share=(a+b)/total
        truth=share[signs].mean()-share[~signs].mean()
        # Units per 10,000 molecules keep the target contrast identical across depths.
        err=w[ok]*share[sel[ok]]
        var=share*(1-share)*(total-k)/(np.maximum(1,total-1)*k)
        mse=(err.sum(1)-truth)**2+(w[ok]**2*var[sel[ok]]).sum(1)
        misses=[hypergeom.pmf(0,total,x,k)[sel].prod(1) for x in (a,b,a+b)]
        q=np.clip(1-misses[0]-misses[1]+misses[2],0,1)
        rmse=10000*np.sqrt(mse.mean())
        mcse=10000*mse.std(ddof=1)/(2*np.sqrt(len(mse)*mse.mean())) if mse.mean()>0 else 0.
        result.append(dict(pair=group,genes=g+'+'+h,contrast_per10k=10000*truth,rmse_per10k=rmse,rmse_mcse=mcse,
                           dual_probability=q.mean(),dual_mcse=q.std(ddof=1)/np.sqrt(len(q)),estimable=ok.sum(),attempts=len(ok)))
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=BASE/'synthetic_pilot.csv')
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--draws',type=int,default=2000)
    ap.add_argument('--section',help='optional single-section demonstration')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    all_data=pd.read_csv(args.input); all_data=all_data[all_data.spatial_cell>=0]
    candidates=[g for g in GENES if g not in ('Acan','Col2a1')]
    abundance=(all_data[candidates].sum()/all_data.total.sum()).sort_values()
    pool=abundance[abundance>0].index.tolist()
    mid=(len(pool)-2)//2
    pairs=dict(primary=('Acan','Col2a1'),rare=tuple(pool[:2]),middle=tuple(pool[mid:mid+2]))
    abundance.rename('pooled_fraction').to_csv(args.out/'abundance.csv')
    rows=[]; original=[]; stress=[]
    for i,(section,d) in enumerate(all_data.groupby('section',sort=True)):
        if args.section and section!=args.section:continue
        # Conservative shared eligibility based on all six panel targets, identical for all budgets.
        named=list(dict.fromkeys(g for pair in pairs.values() for g in pair))
        f=d[d.total-d[named].sum(axis=1)>=100].reset_index(drop=True)
        native=benchmark(d,GENES,12,50,args.draws,24001+23*i,['Acan','Col2a1'])
        native.insert(0,'section',section);original.append(native)
        for m,k in PLANS:
            for j,method in enumerate(METHODS):
                ids,w,ok=policies(f,m,args.draws,31000+i*307+m*11+j,method)
                assert (np.diff(np.sort(ids,axis=1),axis=1)>0).all()
                for row in evaluate(f,ids,w,ok,k,pairs):
                    rows.append(dict(section=section,method=method,sites=m,k=k,eligible=len(f),frame=len(d),**row))
        # Erode glare/dark transition in native spot pitches. A fixed contrast within each
        # support compares designs fairly. This sensitivity changes the estimand explicitly.
        from scipy.spatial import cKDTree
        xy=f[['x','y']].to_numpy()
        pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1])
        all_section=pd.read_csv(args.input)
        buffer=all_section[(all_section.section==section)&(all_section.stratum=='buffer')]
        if len(buffer):
            dist=cKDTree(buffer[['x','y']].to_numpy()).query(xy)[0]/pitch
            for radius in (0,1,2):
                keep=f[dist>radius-1e-7].reset_index(drop=True)
                if keep.stratum.nunique()!=2 or keep.stratum.value_counts().min()<6:continue
                for method in ('random','stratified','stratified_dispersed'):
                    ids,w,ok=policies(keep,12,args.draws,15000+i*53+radius,method)
                    result=evaluate(keep,ids,w,ok,50,{'primary':pairs['primary']})[0]
                    stress.append(dict(section=section,pitches=radius,method=method,eligible=len(keep),**result))
        print(section,flush=True)
    if not original:raise ValueError('No matching section')
    pd.concat(original).to_csv(args.out/'original_five.csv',index=False)
    z=pd.DataFrame(rows)
    z.to_csv(args.out/'fixed600.csv',index=False)
    keys=['section','sites','pair']
    z=z.merge(z[z.method=='random'][keys+['rmse_per10k']].rename(columns={'rmse_per10k':'random_rmse'}),on=keys)
    z['gain']=100*(1-z.rmse_per10k/z.random_rmse)
    z.groupby(['pair','sites','method'])[['gain','rmse_per10k','dual_probability']].median().to_csv(args.out/'fixed600_medians.csv')
    # These retrospective optima use the measured frame and can be tied within MC resolution.
    chosen=[]
    for key,sub in z.groupby(['section','pair','sites']):
        for endpoint,score in [('dual_probability',sub.dual_probability),('rmse_per10k',-sub.rmse_per10k)]:
            best=sub.loc[score.idxmax()]
            chosen.append(dict(section=key[0],pair=key[1],sites=key[2],endpoint=endpoint,method=best.method,estimate=best[endpoint]))
    pd.DataFrame(chosen).to_csv(args.out/'design_lookup.csv',index=False)
    pd.DataFrame(stress).to_csv(args.out/'margin_sensitivity.csv',index=False)
if __name__=='__main__':main()
