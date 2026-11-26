"""Finite-frame panel recovery and resources for repeated whole requests."""
from pathlib import Path
from itertools import combinations
from decimal import Decimal,localcontext
from math import comb,ceil,log1p,isfinite,fsum
import argparse
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.special import betaln

def check_counts(depth,targets,m,k):
    L=np.asarray(depth);Y=np.asarray(targets)
    if L.ndim!=1 or Y.ndim!=2 or len(L)!=len(Y) or not 1<=Y.shape[1]<=6 or not isinstance(k,(int,np.integer)) or k<1 or not isinstance(m,(int,np.integer)) or m<1:
        raise ValueError('Expected a depth vector, one to six target columns, and positive integer m and K.')
    for array in (L,Y):
        if not np.issubdtype(array.dtype,np.number) or not np.isfinite(array).all() or np.any(array<0) or np.any(array!=np.floor(array)):
            raise ValueError('Counts must be finite nonnegative integers.')
    L=L.astype(np.int64);Y=Y.astype(np.int64)
    if np.any(Y.sum(axis=1)>L):raise ValueError('Disjoint target counts exceed total depth.')
    return L,Y

def fixed_request_probability(depth,targets,k,strict=True):
    """All targets recovered in one specified layout, conditional on its selected sites.

    This fixed-layout function gates on total captured molecules. The random-frame
    analysis in evaluate() instead uses the manuscript's panel-excluded depth.
    Short sites contribute their whole library only when strict is false.
    """
    L,Y=check_counts(depth,targets,1,k)
    if len(L)==0 or (strict and np.any(L<k)) or np.any(Y.sum(0)==0):
        return 0.0
    draws=np.minimum(L,k)
    with localcontext() as ctx:
        ctx.prec=65
        result=Decimal(1)
        for size in range(1,Y.shape[1]+1):
            for indices in combinations(range(Y.shape[1]),size):
                miss=Decimal(1)
                for l,y,n in zip(L,Y[:,indices].sum(1),draws):
                    l,y,n=int(l),int(y),int(n)
                    if y>l-n:
                        miss=Decimal(0)
                        break
                    small,large=sorted([y,n])
                    miss*=Decimal(comb(l-large,small))/Decimal(comb(l,small))
                result+=(-1)**size*miss
        return float(max(Decimal(0),min(Decimal(1),result)))

def expected_product(q,m):
    # e_m(q) / C(n,m) from power sums scaled by n**j; m is small.
    n=len(q)
    if n==m:return float(np.prod(q))
    if m>16 or n<2*m:
        # A convex recurrence avoids cancellation in small frames or large requests.
        means=np.zeros(m+1);means[0]=1
        for count,value in enumerate(q,1):
            high=min(m,count);j=np.arange(1,high+1)
            means[1:high+1]=(1-j/count)*means[1:high+1]+(j/count)*value*means[:high]
        return means[m]
    power=q.copy();sums=[];e=[1.]
    for j in range(1,m+1):
        sums.append(float(power.mean())/n**(j-1))
        e.append(fsum((-1)**i*sums[i]*e[-i-1] for i in range(j))/j)
        power*=q
    return e[-1]/(comb(n,m)/n**m)

def decimal_cond(L,Y,k,m):
    # Only used when inclusion-exclusion nearly cancels. Compute from integer
    # frequency classes with extra precision, including zero-depth units.
    with localcontext() as ctx:
        ctx.prec=60
        total=Decimal(1);den=Decimal(comb(len(L),m))
        for s in range(1,Y.shape[1]+1):
            for sub in combinations(range(Y.shape[1]),s):
                triples,freq=np.unique(np.c_[L,Y[:,sub].sum(1),np.minimum(L,k)],axis=0,return_counts=True)
                ps=[Decimal(0)]*m
                for (l,y,n),w in zip(triples,freq):
                    l,y,n=map(int,(l,y,n));n0=min(y,n);other=max(y,n)
                    v=Decimal(comb(l-other,n0))/Decimal(comb(l,n0)) if n0<=l-other else Decimal(0)
                    powv=Decimal(int(w))
                    for j in range(m):
                        powv*=v;ps[j]+=powv
                e=[Decimal(1)]
                for j in range(1,m+1):e.append(sum((-1)**i*ps[i]*e[-i-1] for i in range(j))/j)
                total+=(-1)**s*e[-1]/den
        return float(max(Decimal(0),min(Decimal(1),total)))

def evaluate(depth,targets,m,k,gate=None,precise=True):
    L,Y=check_counts(depth,targets,m,k);N=len(L)
    cutoff=k if gate is None else gate
    if not isinstance(cutoff,(int,np.integer)) or not 1<=cutoff<=k:raise ValueError('gate must be an integer in [1,K]')
    eligible=L-Y.sum(1)>=cutoff;n=int(eligible.sum())
    if n<m:return 0.,np.nan,0.
    att=comb(n,m)/comb(N,m);L=L[eligible];Y=Y[eligible]
    if np.any(Y.sum(0)==0):return att,0.,0.
    draws=np.minimum(L,k);terms=[1.]
    for s in range(1,Y.shape[1]+1):
        for sub in combinations(range(Y.shape[1]),s):
            y=Y[:,sub].sum(1);present=y>0;q=np.ones(n)
            q[present]=hypergeom.pmf(0,L[present],y[present],draws[present])
            terms.append((-1)**s*expected_product(q,m))
    value=fsum(terms)
    if precise and value<1e-8:value=decimal_cond(L,Y,k,m)
    if not -1e-9<=value<=1+1e-9:raise ArithmeticError('Inclusion-exclusion accuracy failure')
    h=min(1.,max(0.,value));return att,h,att*h

def k_threshold(L,Y,m=12,target=.8):
    L,Y=check_counts(L,Y,m,1)
    if not 0<target<1:raise ValueError('Threshold must be in (0,1).')
    N=len(L)
    if N<m:return None,'too few units'
    den=comb(N,m)
    for y in Y.T:
        ceiling=1-comb(int((y==0).sum()),m)/den
        if ceiling<target-1e-12:return None,'native target support below threshold'
    E=L-Y.sum(1);lo=0;hi=int(E.max())+1
    while lo+1<hi:
        mid=(lo+hi)//2
        if comb(int((E>=mid).sum()),m)/den>=target:lo=mid
        else:hi=mid
    def first(left,right):
        if left>right:return None
        # For any K in this interval, eligibility is a subset of eligibility
        # at its left endpoint; recovery under maximal thinning depth dominates.
        upper=evaluate(L,Y,m,right,gate=left,precise=False)[2]
        if upper<target-1e-11:return None
        if left==right:
            p=evaluate(L,Y,m,left)[2]
            if abs(p-target)<1e-10:
                valid=E>=left;cond=decimal_cond(L[valid],Y[valid],left,m)
                p=comb(int(valid.sum()),m)/den*cond
            return left if p>=target else None
        mid=(left+right)//2
        hit=first(left,mid)
        return hit if hit is not None else first(mid+1,right)
    k=first(1,lo)
    return k,'attained' if k is not None else 'no qualifying integer'

def repeat_count(p,t=.8,rho=0,limit=10**12):
    if not (isfinite(p) and 0<=p<=1 and 0<t<1 and 0<=rho<1):raise ValueError('Invalid probability or correlation.')
    if p==0:return None,'no finite count'
    if p>=t:return 1,'finite'
    if rho==0:
        m=ceil(log1p(-t)/log1p(-p));return (m,'finite') if m<=limit else (None,'above numerical limit')
    n=(1-rho)/rho;a=p*n;b=(1-p)*n;goal=log1p(-t)
    logfail=lambda m:betaln(a,b+m)-betaln(a,b)
    lo=0;hi=1
    while hi<limit and logfail(hi)>goal:lo,hi=hi,min(2*hi,limit)
    if logfail(hi)>goal:return None,'above numerical limit'
    while hi-lo>1:
        mid=(lo+hi)//2
        if logfail(mid)<=goal:hi=mid
        else:lo=mid
    return hi,'finite'

def resources(p,budget,area=None,units=None):
    if budget<=0 or (area is not None and area<=0):raise ValueError('Positive budget and area required.')
    o={'p':p,'request_UMIs':budget,'request_mm2':area,'units_per_request':units,
       'expected_independent_requests_to_first':1/p if p else np.nan,
       'expected_nominal_UMIs_to_first':budget/p if p else np.nan,
       'expected_mm2_to_first':area/p if area is not None and p else np.nan}
    for t,rho in [(.8,0),(.95,0),(.8,.1),(.95,.1)]:
        m,status=repeat_count(p,t,rho);key=f'm{int(t*100)}_'+('iid' if rho==0 else 'icc10')
        o[key]=m;o[key+'_status']=status
        o[key+'_planned_UMIs']=m*budget if m is not None else np.nan
        o[key+'_planned_mm2']=m*area if m is not None and area is not None else np.nan
        o[key+'_unit_count']=m*units if m is not None and units is not None else np.nan
    return o

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input_csv',type=Path);ap.add_argument('output_csv',type=Path)
    ap.add_argument('--targets',nargs='+',required=True);ap.add_argument('--depth-col',default='total');ap.add_argument('--sites',type=int,default=12);ap.add_argument('--K',type=int,required=True)
    a=ap.parse_args();df=pd.read_csv(a.input_csv)
    if not 2<=len(a.targets)<=6 or len(a.targets)!=len(set(a.targets)):ap.error('Two to six distinct targets are required.')
    out=[];L=df[a.depth_col].to_numpy()
    for M in range(2,len(a.targets)+1):
        Y=df[a.targets[:M]].to_numpy();att,h,p=evaluate(L,Y,a.sites,a.K);k,status=k_threshold(L,Y,a.sites)
        out.append(dict(M=M,targets='+'.join(a.targets[:M]),K=a.K,A=att,H=h,S=p,K80=k,status=status))
    pd.DataFrame(out).to_csv(a.output_csv,index=False)
