"""Score full requests on a predeclared nominal tissue frame. Spatial coordinates must be registered.
Run beside design_atlas.py from spatial-sampling-budget v0.5.0. All calculations are local.
"""
import argparse,math
from pathlib import Path
import numpy as np,pandas as pd
from design_atlas import METHODS,policies
from scipy.stats import hypergeom

def plan(frame,gene_a,gene_b,budget,site_counts,area,diameter,draws=2000,seed=42,um_per_unit=1.,shape='circle'):
    if not np.isfinite([budget,area,diameter,draws,um_per_unit]).all() or min(budget,area,diameter,draws,um_per_unit)<=0 or gene_a==gene_b:raise ValueError('Set positive resources and two distinct targets.')
    if shape not in ['circle','square']:raise ValueError('Specify circle diameter or an axis-aligned square side.')
    if set(frame.stratum)!= {'dark','light'} or frame.site.duplicated().any():raise ValueError('Supply unique sites in both predefined optical strata.')
    if not set(frame.geometry_ok.dropna().unique())<={0,1} or frame.geometry_ok.isna().any():raise ValueError('geometry_ok must state footprint clearance (0/1) for every site.')
    counts=frame[['total',gene_a,gene_b]].to_numpy(float)
    if not np.isfinite(counts).all() or (counts!=np.floor(counts)).any() or (counts<0).any() or (counts[:,1:].sum(1)>counts[:,0]).any():raise ValueError('Measured counts must be valid integers.')
    L,a,b=counts.astype(np.int64).T;dark=frame.stratum.eq('dark').to_numpy()
    xy=frame[['x','y']].to_numpy(float)*um_per_unit
    if not np.isfinite(xy).all():raise ValueError('Coordinates must be finite.')
    geo=frame.geometry_ok.to_numpy(bool)
    if not (geo&dark&(L>0)).any() or not (geo&~dark&(L>0)).any():raise ValueError('Both strata need at least one measured, in-bound site.')
    freq=np.divide(a+b,L,out=np.zeros(len(L)),where=L>0)
    truth=freq[geo&dark&(L>0)].mean()-freq[geo&~dark&(L>0)].mean()
    rows=[];points=[]
    for m in sorted(set(site_counts)):
        if m<2 or m%2 or budget%m:raise ValueError('Site choices must be even divisors of the budget.')
        k=budget//m
        footprint=m*(math.pi*(diameter/2)**2 if shape=='circle' else diameter**2)/1e6
        for j,method in enumerate(METHODS):
            row=dict(method=method,sites=m,K=k,total_umi=budget,footprint_mm2=footprint,draws=draws)
            if footprint>area or m>len(frame) or (method in ('stratified','stratified_dispersed') and min(dark.sum(),(~dark).sum())<m//2):
                rows.append(dict(row,status='area or nominal capacity exceeded',placement=0.,attainment=0.,strict_dual=0.,failure_risk=1.,estimable=0));continue
            draw_seed=seed+13*m+j
            constructed=np.ones(draws,dtype=bool)
            try:sel,w,both=policies(frame,m,draws,draw_seed,method)
            except ValueError as ex:
                if method!='inhibition':raise
                # Preserve failed proposals in the denominator when some random orders cannot fit.
                sel=np.zeros((draws,m),int);w=np.zeros((draws,m));both=np.zeros(draws,bool)
                seeds=np.random.SeedSequence(draw_seed).generate_state(draws)
                for trial,one_seed in enumerate(seeds):
                    try:
                        ids,weights,valid=policies(frame,m,1,int(one_seed),method)
                        sel[trial],w[trial],both[trial]=ids[0],weights[0],valid[0]
                    except ValueError:
                        constructed[trial]=False
            fit=constructed & geo[sel].all(1)
            for i in range(m-1):
                dv=np.abs(xy[sel[:,i+1:]]-xy[sel[:,i],None,:])
                near=(np.square(dv).sum(2)<diameter**2*(1-1e-12)) if shape=='circle' else (dv<diameter*(1-1e-12)).all(2)
                fit&=~near.any(1)
            full=fit&(L[sel]-a[sel]-b[sel]>=k).all(1)
            ok=full&both
            positive=np.zeros(len(sel));ix=sel[full]
            if len(ix):
                elig=L-a-b>=k;miss=[]
                for y in (a,b,a+b):
                    q=np.ones(len(L));q[elig]=hypergeom.pmf(0,L[elig],y[elig],k)
                    miss.append(q[ix].prod(1))
                positive[full]=np.clip(1-miss[0]-miss[1]+miss[2],0,1)
            if ok.any():
                ix=sel[ok];weights=w[ok]
                v=freq*(1-freq)*np.maximum(L-k,0)/(k*np.maximum(L-1,1))
                error=((weights*freq[ix]).sum(1)-truth)**2+(weights**2*v[ix]).sum(1)
                rmse=10000*np.sqrt(error.mean())
            else:rmse=np.nan
            prob=positive.mean()
            rows.append(dict(row,status='scored',placement=fit.mean(),attainment=full.mean(),full_given_placement=full.sum()/fit.sum() if fit.any() else np.nan,dual_given_full=positive.sum()/full.sum() if full.any() else np.nan,strict_dual=prob,failure_risk=1-prob,estimable=ok.sum(),contrast_RMSE_per10k=rmse,dual_MCSE=positive.std(ddof=1)/np.sqrt(draws) if draws>1 else np.nan))
            if fit.any():
                for site in sel[np.flatnonzero(fit)[0]]:
                    t=frame.iloc[site];points.append(dict(method=method,sites=m,K=k,site=t.site,x=t.x,y=t.y,stratum=t.stratum))
    g=pd.DataFrame(rows)
    if 'contrast_RMSE_per10k' not in g:g['contrast_RMSE_per10k']=np.nan
    mask=(g.placement>0)&g.contrast_RMSE_per10k.notna()
    g['nondominated']=False;g['lowest_error']=False;g['highest_detection']=False
    vals=g.loc[mask,['contrast_RMSE_per10k','strict_dual','footprint_mm2']].to_numpy(copy=True);vals[:,1]*=-1
    if len(vals):
        g.loc[mask,'nondominated']=[not np.any(np.all(vals<=v,axis=1)&np.any(vals<v,axis=1)) for v in vals]
        g.loc[mask,'lowest_error']=g.loc[mask,'contrast_RMSE_per10k']==vals[:,0].min()
    if (g.placement>0).any() and g.strict_dual.max()>0:g['highest_detection']=(g.placement>0)&g.strict_dual.eq(g.strict_dual.max())
    return g,pd.DataFrame(points,columns=['method','sites','K','site','x','y','stratum'])
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--section',required=True);ap.add_argument('--targets',nargs=2,required=True);ap.add_argument('--budget',type=int,required=True);ap.add_argument('--sites',type=int,nargs='+',required=True);ap.add_argument('--area-mm2',type=float,required=True);ap.add_argument('--size-um',type=float,required=True);ap.add_argument('--shape',choices=['circle','square'],default='circle');ap.add_argument('--um-per-unit',type=float,required=True);ap.add_argument('--draws',type=int,default=2000);ap.add_argument('--seed',type=int,default=42);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();f=pd.read_csv(a.input);f=f[(f.section==a.section)&f.stratum.isin(['dark','light'])].reset_index(drop=True)
    g,loc=plan(f,*a.targets,a.budget,a.sites,a.area_mm2,a.size_um,a.draws,a.seed,a.um_per_unit,a.shape)
    a.out.mkdir(parents=True,exist_ok=True);g.to_csv(a.out/'protocols.csv',index=False);loc.to_csv(a.out/'first_feasible_sites.csv',index=False)
if __name__=='__main__':main()
