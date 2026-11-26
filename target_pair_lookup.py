"""Find the first K with 80% strict dual recovery in a 12-unit request.
Test all integers that can satisfy the budget-attainment bound; exact without-replacement scoring.
"""
from pathlib import Path
import math,json
import numpy as np,pandas as pd
from scipy.stats import hypergeom
from extraction_bias_sim import repeat_count

pairs={'Matrix synthesis':('ACAN','COL2A1'),'Matrix turnover':('MMP3','ADAMTS5'),'Collagenase / cytokine':('MMP13','IL1B'),'Notochordal markers':('KRT19','TBXT')}
M=12
def sym(v):
    n=len(v)
    if n<M:return 0.
    # Newton identities on normalized power sums. The terminal scale retains the
    # exact without-replacement denominator for a uniform distinct-unit request.
    v=np.asarray(v,float);power=[];prod=np.ones_like(v)
    for j in range(1,M+1):
        prod*=v;power.append(prod.mean()/n**(j-1))
    seq=[1.]
    for j in range(1,M+1):
        seq.append(sum((-1)**i*power[i]*seq[j-1-i] for i in range(j))/j)
    return seq[-1]/(math.comb(n,M)/n**M)
def passp(n,N):
    return math.comb(n,M)/math.comb(N,M) if n>=M else 0.
def exact(L,a,b,k):
    ok=L-a-b>=k;n=ok.sum();P=passp(n,len(L))
    if not P:return 0.,0.,'no full request'
    if not a[ok].any() or not b[ok].any():return P,0.,'target absent among eligible units'
    q=[hypergeom.pmf(0,L[ok],g[ok],k) for g in (a,b,a+b)]
    dual=1-sym(q[0])-sym(q[1])+sym(q[2])
    # Small positive recovery must not be mislabeled an impossible event.
    if dual<1e-9:
        states=np.array([q[2],q[1]-q[2],q[0]-q[2],np.maximum(0,1-q[0]-q[1]+q[2])]).T
        dp=np.zeros((4,M+1));dp[0,0]=1.
        for i,(no,x,y,xy) in enumerate(states,1):
            js=np.arange(1,min(i,M)+1);old=dp[:,js-1]
            fresh=np.stack([no*old[0],no*old[1]+x*(old[0]+old[1]),no*old[2]+y*(old[0]+old[2]),old[3]+x*old[2]+y*old[1]+xy*(old[0]+old[1]+old[2])])
            dp[:,js]=(1-js/i)*dp[:,js]+js/i*fresh
        dual=dp[-1,-1]
    return P,float(P*np.clip(dual,0,1)),'evaluated'
def one(z,kind,unit,label,A,B,depth):
    row=dict(source=kind,unit=unit,group=label,target_a=A,target_b=B,sites=M,reference_K=depth,frame=len(z))
    if z[[A,B]].isna().any().any():return dict(row,status='feature not represented')
    L,a,b=z[['total',A,B]].to_numpy(np.int64).T;N=len(L);e=L-a-b
    if not (e>=0).all() or N<M: raise ValueError('Provide at least 12 distinct sites with compatible totals.')
    lo,hi=0,int(e.max())+1
    while lo+1<hi:
        k=(lo+hi)//2
        if passp(sum(e>=k),N)>=.8:lo=k
        else:hi=k
    # Even without thinning, a uniform request on this fixed frame can
    # recover both more often than this full-depth pooled detection probability.
    native=1-passp(np.sum(a==0),N)-passp(np.sum(b==0),N)+passp(np.sum((a+b)==0),N)
    pk,sk,msg=exact(L,a,b,depth)
    m,ms=repeat_count(sk,.8,0);mc,mcstatus=repeat_count(sk,.8,.10)
    k80=np.nan;p_at=np.nan;p_prev=np.nan;ntry=0
    if native<.8-1e-12:why='native target support below 80%'
    else:
        why='no qualifying integer before attainment falls below 80%'
        for k in range(1,lo+1):
            _,s,_=exact(L,a,b,k);ntry+=1
            if s>=.8:
                k80=k;p_at=s;why='attained';break
            p_prev=s
    return dict(row,native_dual_ceiling=max(0,native),max_K_with_80pct_attainment=lo,integers_scored=ntry,K80_per_unit=k80,K80_request_total=12*k80 if np.isfinite(k80) else np.nan,
        status=why,p_at_K80=p_at,p_previous_to_K80=p_prev if np.isfinite(k80) else np.nan,
        full_at_reference=pk,strict_at_reference=sk,reference_status=msg,m80_at_reference=m,repeat_status=ms,
        m80_rho010=mc,correlated_repeat_status=mcstatus,abundance_A=a.sum()/L.sum(),abundance_B=b.sum()/L.sum())
def main():
    import argparse
    ap=argparse.ArgumentParser(description='Uniform requests of 12 distinct sites with a fixed per-site non-target gate.')
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--unit-column',default='section')
    ap.add_argument('--targets',nargs=2,required=True)
    ap.add_argument('--reference-K',type=int,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ar=ap.parse_args();a,b=ar.targets
    if a==b or ar.reference_K<1:ap.error('Specify two distinct targets and positive K.')
    d=pd.read_csv(ar.input)
    if ar.unit_column not in d or d[ar.unit_column].isna().any():ap.error('Missing unit identifiers.')
    if not {a,b,'total','site'}.issubset(d) or d.site.isna().any() or d.duplicated([ar.unit_column,'site']).any():ap.error('Supply unique site identifiers and count columns.')
    x=d[[a,b,'total']].to_numpy(float)
    if not np.isfinite(x).all() or (x<0).any() or (x!=np.floor(x)).any() or (x[:,0]+x[:,1]>x[:,2]).any():ap.error('Counts must be nonnegative integers with a compatible total.')
    rows=[one(z,'supplied frame',str(k),'user targets',a,b,ar.reference_K) for k,z in d.groupby(ar.unit_column)]
    pd.DataFrame(rows).to_csv(ar.out,index=False)
if __name__=='__main__':main()
