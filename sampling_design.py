"""Portable planning utilities for budget-limited tissue assays."""
import argparse
import json
from math import ceil, log1p, comb
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import betaln
from scipy.spatial.distance import cdist

POLICIES=("random","systematic_grid","stratified","stratified_dispersed","spatial_only")


def repeat_plan(p,target,rho=0):
    if not 0<=p<=1 or not 0<target<1 or not 0<=rho<1:
        raise ValueError("require 0 ≤ p ≤ 1, 0 < target < 1 and 0 ≤ rho < 1")
    if p==0:return None
    if p==1:return 1
    if rho==0:return ceil(log1p(-target)/log1p(-p))
    a=p*(1/rho-1); b=(1-p)*(1/rho-1)
    def log_failure(n):return betaln(a,b+n)-betaln(a,b)
    lo,hi=0,1
    while log_failure(hi)>log1p(-target):
        if hi>=2**40:raise ValueError("required count exceeds supported numerical search range")
        hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if log_failure(mid)<=log1p(-target):hi=mid
        else:lo=mid
    return hi


def recover(total,a,b,k):
    if min(a,b,total,k)<0 or a+b>total or k<1:raise ValueError("invalid molecule counts or budget")
    if a==0 or b==0:return 0.
    n=min(total,k)
    denominator=comb(total,n)
    def miss(g):return comb(total-g,n)/denominator if total-g>=n else 0.
    return max(0.,min(1.,1-miss(a)-miss(b)+miss(a+b)))


def compare_blocks(atoms,blocks,genes,k):
    if genes[0]==genes[1] or k<1:raise ValueError("two distinct targets and a positive integer budget required")
    if atoms.site.duplicated().any() or blocks.block.duplicated().any():raise ValueError("duplicate ID")
    counts=atoms.set_index("site")[["total",*genes]]
    vals=counts.to_numpy(dtype=float)
    if not np.isfinite(vals).all() or (vals<0).any() or (vals!=np.floor(vals)).any() or (vals[:,1:].sum(1)>vals[:,0]).any():
        raise ValueError("atomic counts must be finite, nonnegative integers with targets included in total")
    sets=[str(ids).split("|") for ids in blocks.sites]
    if len(set(map(len,sets)))!=1 or any(len(set(s))!=len(s) for s in sets):raise ValueError("use equal-sized complete blocks")
    rows=[]
    for i,j in combinations(range(len(blocks)),2):
        arm=" + ".join(sorted([blocks.zone.iloc[i],blocks.zone.iloc[j]]))
        if set(sets[i]).intersection(sets[j]):
            rows.append(dict(block1=blocks.block.iloc[i],block2=blocks.block.iloc[j],layout=arm,fit=0,full=0,strict=0.,short_accepted=0.))
            continue
        vals=counts.loc[sets[i]+sets[j]].sum().to_numpy(dtype=float)
        if np.any(vals!=np.floor(vals)):raise ValueError("integer counts required")
        total,a,b=map(int,vals)
        q=recover(total,a,b,k)
        rows.append(dict(block1=blocks.block.iloc[i],block2=blocks.block.iloc[j],
                         layout=arm,fit=1,
                         full=int(total>=k),strict=q if total>=k else 0.,short_accepted=q))
    return pd.DataFrame(rows,columns=["block1","block2","layout","fit","full","strict","short_accepted"])


def calibration(table):
    required=["observed","source_expected","per_library_multiplier","common_multiplier"]
    vals=table[required].to_numpy(dtype=float)
    if not np.isfinite(vals).all() or (vals<0).any() or (vals[:,1]<=0).any():raise ValueError("invalid scale inputs")
    if table.common_multiplier.nunique()!=1 or table.duplicated(["specimen","gene"]).any():raise ValueError("check calibration/row keys")
    if table.groupby("specimen").per_library_multiplier.nunique().max()!=1:raise ValueError("per-library factor must be fixed within specimen")
    for name,factor in [("unscaled",1),("per_library",table.per_library_multiplier),("common",table.common_multiplier)]:
        table[name]=abs(table.observed-table.source_expected*factor)/table.source_expected
    return table.groupby("specimen")[["unscaled","per_library","common"]].mean()


def partition(indices,xy,k):
    if len(indices)<k:raise ValueError("insufficient eligible sites per stratum")
    if k==1:return [indices]
    axis=np.argmax(np.ptp(xy[indices],axis=0))
    ids=indices[np.argsort(xy[indices,axis],kind="stable")]
    a=k//2; n=int(round(len(ids)*a/k))
    return partition(ids[:n],xy,a)+partition(ids[n:],xy,k-a)


def draws(frame,m,n,seed,method):
    if method not in POLICIES or m<2 or m%2 or n<1 or len(frame)<m:raise ValueError("invalid design")
    if frame.site.duplicated().any():raise ValueError("duplicate site")
    rng=np.random.default_rng(seed)
    xy=frame[["x","y"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():raise ValueError("non-finite coordinates")
    strata=frame.stratum.to_numpy()
    if set(strata)!={"light","dark"}:raise ValueError("two image-defined strata required")
    idx=np.arange(len(frame)); sign=np.where(strata=="dark",1.,-1.)
    pop=[idx[sign==s] for s in [1,-1]]
    if method in ("stratified","stratified_dispersed") and min(map(len,pop))<m//2:
        raise ValueError("insufficient eligible sites per stratum")
    if method=="random":
        chosen=np.stack([rng.choice(idx,m,replace=False) for _ in range(n)])
    elif method=="stratified":
        chosen=np.stack([np.concatenate([rng.choice(ids,m//2,replace=False) for ids in pop]) for _ in range(n)])
    elif method in ("stratified_dispersed","spatial_only"):
        if method=="stratified_dispersed" and "spatial_cell" in frame and all(frame.iloc[g].spatial_cell.nunique()==m//2 for g in pop):
            bins=[g[frame.iloc[g].spatial_cell.to_numpy()==v] for g in pop for v in sorted(frame.iloc[g].spatial_cell.unique())]
        elif method=="stratified_dispersed":
            bins=[cell for g in pop for cell in partition(g,xy,m//2)]
        else:bins=partition(idx,xy,m)
        chosen=np.column_stack([rng.choice(g,n) for g in bins])
    else:
        # Randomly translated rectangular lattice in the bounding box. Snap each
        # lattice node to the nearest unused tissue site in random node order.
        factors=[(i,m//i) for i in range(1,m+1) if m%i==0]
        span=np.ptp(xy,axis=0)
        nx,ny=min(factors,key=lambda q:abs(np.log((span[0]+1e-9)/(span[1]+1e-9)*q[1]/q[0])))
        cell=span/[nx,ny]
        nodes=np.stack(np.meshgrid(np.arange(nx),np.arange(ny)),axis=-1).reshape(-1,2)
        chosen=np.empty((n,m),dtype=int)
        for j in range(n):
            positions=xy.min(0)+(nodes+rng.random(2))*cell
            dist=cdist(positions[rng.permutation(m)],xy)
            for z in range(m):
                loc=np.argmin(dist[z]);chosen[j,z]=loc;dist[:,loc]=np.inf
    if method in ("stratified_dispersed","spatial_only"):base=np.array([len(i) for i in bins],float)[None,:]
    else:base=np.ones((1,m))
    a=(sign[chosen]>0)*base;b=(sign[chosen]<0)*base
    ok=(a.sum(1)>0)&(b.sum(1)>0)
    w=a/np.maximum(1,a.sum(1))[:,None]-b/np.maximum(1,b.sum(1))[:,None]
    return chosen,w,ok


def benchmark(frame,genes,m,k,n,seed,eligibility_targets=None):
    if k<=0 or n<=0 or len(genes)!=len(set(genes)):raise ValueError("invalid budget, draw count or genes")
    y=frame[genes].to_numpy(float)
    total=frame.total.to_numpy(float)
    if not np.isfinite(y).all() or not np.isfinite(total).all() or np.any(y<0) or np.any(total<0) or np.any(y.sum(1)>total) or np.any(y!=np.floor(y)) or np.any(total!=np.floor(total)):
        raise ValueError("counts must be finite nonnegative integers")
    primary=eligibility_targets or genes
    if len(set(primary))!=len(primary) or not set(primary)<=set(genes):raise ValueError("eligibility targets must be unique measured genes")
    target_total=frame[primary].sum(axis=1).to_numpy()
    eligible=total-target_total>=k
    labels=list(genes)
    if len(primary)>1:
        y=np.column_stack((y,target_total))
        labels.append("+".join(primary))
    f=frame.loc[eligible].reset_index(drop=True)
    L=total[eligible,None]; fraction=y[eligible]/L
    expected=k*fraction
    variance=k*fraction*(1-fraction)*(L-k)/np.maximum(1,L-1)
    truth=expected[f.stratum=="dark"].mean(0)-expected[f.stratum=="light"].mean(0)
    rows=[]
    for j,method in enumerate(POLICIES):
        ids,w,ok=draws(f,m,n,seed+j,method)
        if not ok.any():
            rows.extend(dict(method=method,gene=g,frame=len(frame),eligible=len(f),estimable=0,draws=n,rmse=None) for g in labels)
            continue
        w=w[ok];ids=ids[ok]
        residual=np.einsum("ij,ijk->ik",w,expected[ids])-truth
        mse=residual**2+np.einsum("ij,ijk->ik",w**2,variance[ids])
        for gene,r in zip(labels,np.sqrt(mse.mean(0))):
            rows.append(dict(method=method,gene=gene,frame=len(frame),eligible=len(f),estimable=int(ok.sum()),draws=n,rmse=r))
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__); sub=ap.add_subparsers(dest="action",required=True)
    a=sub.add_parser("budget");a.add_argument("--p",type=float,required=True);a.add_argument("--rho",type=float,default=0)
    a=sub.add_parser("sites");a.add_argument("frame");a.add_argument("--method",choices=(*POLICIES,"maximin","inhibition"),default="stratified_dispersed")
    a.add_argument("--m",type=int,default=12);a.add_argument("--seed",type=int,default=42)
    a=sub.add_parser("calibrate");a.add_argument("table")
    a=sub.add_parser("layouts");a.add_argument("atoms");a.add_argument("blocks");a.add_argument("--targets",nargs=2,required=True);a.add_argument("--k",type=int,required=True)
    a.add_argument("--credit",type=float,default=0.)
    a=sub.add_parser("benchmark");a.add_argument("frame");a.add_argument("--genes",nargs="+",required=True)
    a.add_argument("--eligibility-targets",nargs="+")
    a.add_argument("--k",type=int,required=True);a.add_argument("--m",type=int,default=12);a.add_argument("--draws",type=int,default=1000);a.add_argument("--seed",type=int,default=42)
    a=ap.parse_args()
    if a.action=="budget":print(json.dumps({f"m{round(t*100)}":repeat_plan(a.p,t,a.rho) for t in [.8,.95]}))
    elif a.action=="sites":
        f=pd.read_csv(a.frame,dtype={"site":str})
        from design_atlas import policies
        ix,_,ok=policies(f,a.m,1,a.seed,a.method)
        print(f.iloc[ix[0]][["site","x","y","stratum"]].assign(contrast_estimable=bool(ok[0])).to_csv(index=False),end="")
    elif a.action=="calibrate":print(calibration(pd.read_csv(a.table)).to_csv(),end="")
    elif a.action=="benchmark":
        print(benchmark(pd.read_csv(a.frame,dtype={"site":str}),a.genes,a.m,a.k,a.draws,a.seed,a.eligibility_targets).to_csv(index=False),end="")
    else:
        if not 0<=a.credit<=1:raise ValueError("credit must lie in [0,1]")
        r=compare_blocks(pd.read_csv(a.atoms,dtype={"site":str}),pd.read_csv(a.blocks,dtype={"sites":str}),a.targets,a.k)
        summary=r.groupby("layout")[["fit","full","strict","short_accepted"]].mean()
        summary["planned_pairs"]=r.groupby("layout").size()
        summary["legal_pairs"]=r.groupby("layout").fit.sum()
        summary["strict_given_placement"]=summary.strict.div(summary.fit.replace(0,np.nan))
        summary["dual_given_full"]=summary.strict.div(summary.full.replace(0,np.nan))
        summary["utility"]=summary.strict+a.credit*(summary.short_accepted-summary.strict)
        summary["highest_utility"]=np.isclose(summary.utility,summary.utility.max(),atol=1e-12,rtol=0)
        print(summary.to_csv(),end="")


if __name__=="__main__":main()
