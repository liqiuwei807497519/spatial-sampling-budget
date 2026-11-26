"""Regional extraction sensitivity for two disjoint targets in registered rectangular blocks.

Empirical input columns: reader, section (capture), profile, zone (central/lamellar),
row, col, height, width, L, a, c. All coordinates are on a common integer bin lattice.
The primary result integrates binomial retention. A seeded 200-realisation estimator
independently samples the same retention distribution and reuses each pool over all K.
"""
import argparse,math
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.stats import binom,nbinom
from scipy.special import betaln

KS=np.arange(50,401);RATIOS=[.4,.5,.6,.7,.8,.9,1.,1.15,1.3,1.5,1.75,2.]
LIMIT=10**12
def repeat_count(p,target,icc):
    if p<=0:return np.nan,'zero probability'
    if p>=target:return 1,'finite'
    goal=np.log1p(-target)
    if icc==0:
        m=math.ceil(goal/np.log1p(-p));return (m,'finite') if m<=LIMIT else (np.nan,'above 1e12')
    a=p*(1/icc-1);b=(1-p)*(1/icc-1);norm=betaln(a,b)
    def enough(m):return betaln(a,b+m)-norm<=goal
    if not enough(LIMIT):return np.nan,'above 1e12'
    lo,hi=0,1
    while not enough(hi):lo,hi=hi,min(2*hi,LIMIT)
    while lo+1<hi:
        m=(lo+hi)//2
        if enough(m):hi=m
        else:lo=m
    return hi,'finite'
def repeats(d):
    cols={}
    for arm in ['across','lateral','best']:
        values=d['S_'+arm].to_numpy()
        for target in (.8,.95):
            for icc in (0.,.1):
                key=f'm{round(100*target)}_{arm}_'+('ind' if not icc else 'rho010')
                result=[repeat_count(p,target,icc) for p in values]
                cols[key]=[r[0] for r in result];cols[key+'_status']=[r[1] for r in result]
    return pd.concat([d,pd.DataFrame(cols,index=d.index)],axis=1)
def region_eta(r):
    if not np.isfinite(r) or r<=0:raise ValueError('rho_ext must be positive')
    return min(1.,r),min(1.,1/r)
def pool_hist(counts,weights=None):
    L,a,c=counts.T;h=[]
    for y in (a,c,a+c):
        shape=(int(L.max())+1,int(y.max())+1)
        h.append(np.bincount(L*shape[1]+y,weights=weights,minlength=np.prod(shape)).reshape(shape).astype(float))
    z=h[0].sum()
    return [m/z for m in h]
def combine(left,right):
    # Cross-region block libraries are independent draws from disjoint regional pools.
    h=[]
    for x,y in zip(left,right):
        a=fftconvolve(x,y);a[a<2e-17]=0;h.append(a/a.sum())
    return h
def thin_hist(h,eta):
    if eta==1:return h
    @lru_cache(None)
    def pmf(n):return binom.pmf(np.arange(n+1),n,eta)
    out=[]
    for before in h:
        after=np.zeros_like(before)
        for L,y in zip(*np.nonzero(before)):
            py=pmf(y)*before[L,y];pr=pmf(L-y)
            for j,weight in enumerate(py):after[j:j+len(pr),j]+=weight*pr
        np.testing.assert_allclose(after.sum(),1,atol=1e-11);out.append(after)
    return out
def mc_hist(counts,eta,rep,rng,w=None):
    if eta==1:return pool_hist(counts,w)
    L,a,c=counts.T;shape=(rep,len(L))
    A=rng.binomial(a,eta,size=shape);C=rng.binomial(c,eta,size=shape)
    Lnew=A+C+rng.binomial(L-a-c,eta,size=shape)
    weights=np.tile(w,rep) if w is not None else None
    return pool_hist(np.column_stack([Lnew.ravel(),A.ravel(),C.ravel()]),weights)
def joint_by_depth(h):
    # For each K and each original L >= K, aggregate the unconditional dual mass at that L.
    dep=h[0].sum(axis=1);signed=np.zeros((len(dep),max(a.shape[1] for a in h)))
    for sign,mat in zip([1,1,-1],h):signed[:,:mat.shape[1]]+=sign*mat
    idx,y=np.nonzero(np.abs(signed[:,1:])>1e-17);y+=1;weight=signed[idx,y];q=np.ones(len(y));out=np.zeros((401,len(dep)))
    for k in range(1,401):
        den=np.maximum(idx-k+1,1)
        q*=np.maximum(idx-k+1-y,0)/den
        out[k]=np.bincount(idx,weights=weight*(1-q),minlength=len(dep));out[k,:min(k,len(dep))]=0
    return dep,np.maximum(0,out)
def direct_curve(h):
    dep,j=joint_by_depth(h);depth=np.arange(len(dep))
    pos=dep-h[0][:,0]-h[1][:,0]+h[2][:,0]
    u=np.r_[0,np.maximum(pos,0).cumsum()]
    return np.column_stack([np.array([dep[k:].sum() for k in KS]),j[KS].sum(1),u[np.minimum(KS,len(dep))]])
def uniform_retention(dep,j,eta):
    # Given retained depth >=K, a K-subset of uniformly thinned molecules is a uniform
    # K-subset of the original pair. Integrate support and short samples over all depths.
    if eta==1:
        full=(np.arange(len(dep))[None,:]>=KS[:,None]);pm=np.where(np.arange(len(dep))[None,:]==np.arange(400)[:,None],1.,0.)
    else:
        F=np.arange(len(dep))[None,:]
        full=binom.sf(KS[:,None]-1,F,eta)
        pm=binom.pmf(np.arange(400)[:,None],F,eta)
    short=np.r_[0,(pm*j[:400]).sum(1).cumsum()]
    return np.column_stack([full@dep,(j[KS]*full).sum(1),short[KS]])
def legal_totals(blk):
    xy=blk[['row','col']].to_numpy(int);sz=blk[['height','width']].to_numpy(int)
    if np.any(sz<=0):raise ValueError('Invalid rectangle')
    a=blk[['L','a','c']].to_numpy(int)
    i,j=np.triu_indices(len(a),1)
    legal=(np.minimum(xy[i]+sz[i],xy[j]+sz[j])<=np.maximum(xy[i],xy[j])).any(1)
    v=a[i[legal]]+a[j[legal]]
    # An integer radix avoids a slow sort of N x 3 rows.
    radix=int(max(v[:,1:].max()+1,2));key=(v[:,0]*radix+v[:,1])*radix+v[:,2]
    keys,w=np.unique(key,return_counts=True);return np.column_stack([keys//radix//radix,keys//radix%radix,keys%radix]),w
def per_profile(group,ratios,mc,seed):
    center=group.query("zone=='central'");side=group.query("zone=='lamellar'")
    v=[s[['L','a','c']].to_numpy(int) for s in (center,side)]
    if any(np.any(x[:,0]<x[:,1]+x[:,2]) or (x<0).any() for x in v):raise ValueError('Invalid counts')
    cxy=center[['row','col']].to_numpy();sxy=side[['row','col']].to_numpy()
    chi=cxy+center[['height','width']].to_numpy();shi=sxy+side[['height','width']].to_numpy()
    for xy,hi in zip(cxy,chi):
        if ((np.minimum(hi,shi)>np.maximum(xy,sxy)).all(1)).any():raise ValueError('Regions overlap at native bins')
    lat,mass=legal_totals(side);lh=pool_hist(lat,mass);dep,j=joint_by_depth(lh)
    hc,ho=map(pool_hist,v);rows=[]
    for ir,rho in enumerate(ratios):
        eo,ec=region_eta(rho);rng=np.random.default_rng(seed+ir)
        exact=direct_curve(combine(thin_hist(hc,ec),thin_hist(ho,eo)));ll=uniform_retention(dep,j,eo)
        if rho==1:ac_mc,lat_mc=exact,ll
        else:
            ac_mc=direct_curve(combine(mc_hist(v[0],ec,mc,rng),mc_hist(v[1],eo,mc,rng)))
            lat_mc=direct_curve(mc_hist(lat,eo,mc,rng,mass)) if eo<1 else ll
        for estimator,ac,la in [('integrated',exact,ll),('MC200' if mc==200 else 'MC'+str(mc),ac_mc,lat_mc)]:
            p,S,B=ac.T;l,T,D=la.T
            h=np.divide(S,p,out=np.full_like(p,np.nan),where=p>1e-14)
            g=np.divide(T,l,out=np.full_like(l,np.nan),where=l>1e-14)
            delta=100*(S-T);delivery=50*(p-l)*(h+g);conditional=50*(h-g)*(p+l)
            valid=np.isfinite(delivery)
            np.testing.assert_allclose((delivery+conditional)[valid],delta[valid],atol=5e-10)
            rows.append(pd.DataFrame(dict(rho_ext=rho,K=KS,estimator=estimator,p_across=p,p_lateral=l,S_across=S,S_lateral=T,
                B_across=B,B_lateral=D,net_pp=delta,delivery_pp=np.nan_to_num(delivery),conditional_pp=np.nan_to_num(conditional),
                unsupported_pp=np.where(valid,0,delta),supported=valid.astype(int))))
    return pd.concat(rows,ignore_index=True)
def empirical(path,out,mc=200):
    df=pd.read_csv(path);allrows=[]
    for i,(keys,grp) in enumerate(df.groupby(['reader','section','profile'],sort=True)):
        vals=per_profile(grp,RATIOS,mc,20261126+57*i)
        for key,v in zip(('outline','capture','profile'),keys):vals[key]=v
        allrows.append(vals);print('profile',i+1,keys,flush=True)
    profile=pd.concat(allrows,ignore_index=True);profile.to_csv(out/'extraction_profiles.csv.gz',index=False)
    g=profile.groupby(['outline','capture','rho_ext','estimator','K'],as_index=False).mean(numeric_only=True)
    sizes=profile.groupby(['outline','capture']).profile.nunique();assert (sizes==4).all()
    return finish(g,out,'empirical')
def generative(out,p=.003):
    rows=[];k=KS;counts=np.arange(400)
    def both(n,f):return np.maximum(0,1-2*np.exp(n*np.log1p(-f))+np.exp(n*np.log1p(-2*f)))
    for rho in RATIOS:
        eo,ec=region_eta(rho);L=180*eo;C=540*ec
        val=[]
        for mean,f in [(C+L,p*(.45*C+L)/(C+L)),(2*L,p)]:
            prob=18/(18+mean);att=nbinom.sf(k-1,18,prob);h=both(k,f)
            short=np.r_[0,np.cumsum(nbinom.pmf(counts,18,prob)*both(counts,f))][k]
            val.append((att,h,att*h,short))
        (a,h,s,b),(l,g,t,d)=val
        rows.append(pd.DataFrame(dict(outline='model',capture='model',rho_ext=rho,estimator='integrated',K=k,p_across=a,p_lateral=l,S_across=s,S_lateral=t,B_across=b,B_lateral=d,net_pp=100*(s-t),
            delivery_pp=50*(a-l)*(h+g),conditional_pp=50*(h-g)*(a+l),unsupported_pp=0.,supported=1.,outer_fraction=p,shape=18)))
    return finish(pd.concat(rows,ignore_index=True),out,'generative')
def finish(g,out,mode):
    g['eta_outer']=np.minimum(1,g.rho_ext);g['eta_central']=np.minimum(1,1/g.rho_ext)
    g['credit_slope_pp']=100*(g.B_across-g.B_lateral);g['net_full_credit_pp']=g.net_pp+g.credit_slope_pp
    g['S_best']=g[['S_across','S_lateral']].max(axis=1)
    np.testing.assert_allclose(g.net_pp,g.delivery_pp+g.conditional_pp+g.unsupported_pp,atol=3e-10)
    g.to_csv(out/f'extraction_{mode}_grid.csv',index=False)
    selects=g[g.K.isin([173,259])].copy();repeats(selects).to_csv(out/f'extraction_{mode}_repeats.csv',index=False)
    intervals=[]
    for group,z in g.groupby(['outline','capture','rho_ext','estimator']):
        z=z.sort_values('K')
        for w in (0.,1.):
            delta=z.net_pp+w*z.credit_slope_pp
            labels=np.where(delta>1e-10,'across',np.where(delta< -1e-10,'lateral','tie'))
            edges=np.r_[0,np.flatnonzero(labels[1:]!=labels[:-1])+1,len(z)]
            for low,high in zip(edges[:-1],edges[1:]):
                intervals.append(dict(zip(['outline','capture','rho_ext','estimator'],group),w=w,start=int(z.K.iloc[low]),end=int(z.K.iloc[high-1]),preferred=labels[low]))
    pd.DataFrame(intervals).to_csv(out/f'extraction_{mode}_intervals.csv',index=False)
    return g
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mode',choices=['empirical','generative'],required=True)
    ap.add_argument('--blocks',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--n-mc',type=int,default=200)
    args=ap.parse_args();args.out.mkdir(exist_ok=True,parents=True)
    if args.mode=='empirical':
        if args.blocks is None or args.n_mc<1:ap.error('Supply registered blocks and a positive number of MC realisations')
        empirical(args.blocks,args.out,args.n_mc)
    else:generative(args.out)
