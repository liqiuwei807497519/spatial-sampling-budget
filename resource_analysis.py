"""Compare fixed footprints and molecular budgets, with design-based weighting sensitivity."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from design_atlas import METHODS, PLANS, policies, evaluate
from sampling_design import partition

def nondominated(points):
    points=np.asarray(points,float)
    if not np.isfinite(points).all():raise ValueError("All objectives must be finite")
    return np.array([not np.any(np.all(points<=row,axis=1)&np.any(points<row,axis=1)) for row in points])

def ht_weights(frame,ids):
    dark=frame.stratum.eq('dark').to_numpy()
    N1=int(dark.sum());N0=len(frame)-N1
    cells=partition(np.arange(len(frame)),frame[['x','y']].to_numpy(),ids.shape[1])
    inverse_pi=np.zeros(len(frame))
    for j,cell in enumerate(cells):
        if not np.isin(ids[:,j],cell).all():raise ValueError('draw and cell assignments differ')
        inverse_pi[cell]=len(cell)
    return inverse_pi[ids]*np.where(dark[ids],1/N1,-1/N0)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=Path(__file__).parent/'synthetic_pilot.csv')
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--fixed600',type=Path,required=True)
    args=ap.parse_args();args.out.mkdir(exist_ok=True,parents=True)
    master=pd.read_csv(args.fixed600)
    z=master.query("pair=='primary'").copy()
    spot_area=np.pi*(.055/2)**2
    z['footprint_mm2']=z.sites*spot_area
    z['area_vs_six']=z.sites/6
    z['total_umi']=z.sites*z.k
    cols=['rmse_per10k','dual_probability','footprint_mm2']
    def front(g):
        values=g[cols].to_numpy().copy();values[:,1]*=-1
        return pd.Series(nondominated(values),index=g.index)
    z['nondominated']=pd.concat([front(g) for _,g in z.groupby('section')])
    z.to_csv(args.out/'Pareto_sections.csv',index=False)
    med=z.groupby(['method','sites','k','area_vs_six','footprint_mm2'],as_index=False)[['rmse_per10k','dual_probability','nondominated']].agg({'rmse_per10k':'median','dual_probability':'median','nondominated':'sum'}).rename(columns={'nondominated':'sections_on_front'})
    med['median_front']=front(med)
    med.to_csv(args.out/'Pareto_medians.csv',index=False)
    data=pd.read_csv(args.input).query('spatial_cell>=0')
    rows=[];areas=[];pair={'primary':('Acan','Col2a1')}
    names=list(dict.fromkeys(g for pair in master.genes.unique() for g in pair.split('+')))
    for i,(s,section) in enumerate(data.groupby('section',sort=True)):
        frame=section[section.total-section[names].sum(axis=1)>=100].reset_index(drop=True)
        for m,k in PLANS:
            for j,method in enumerate(METHODS):
                if m!=12 and method!='spatial_only':continue
                seed=31000+i*307+m*11+j
                sel,w,ok=policies(frame,m,2000,seed,method)
                def record(weights,mask,depth):
                    return evaluate(frame,sel,weights,mask,depth,pair)[0]
                baseline=record(w,ok,k)
                expected=z[(z.section==s)&(z.method==method)&(z.sites==m)].iloc[0]
                np.testing.assert_allclose([baseline['rmse_per10k'],baseline['dual_probability']],[expected.rmse_per10k,expected.dual_probability],atol=1e-11,rtol=0)
                if method=='spatial_only':
                    ht=ht_weights(frame,sel)
                    for tag,weights,mask in [('Hajek_paired',w,ok),('HT_paired',ht,ok),('HT_all',ht,np.ones(len(ok),bool))]:
                        v=record(weights,mask,k);scores=z[(z.section==s)&(z.sites==m)&(z.method!=method)].rmse_per10k.to_numpy()
                        rows.append(dict(section=s,sites=m,k=k,estimator=tag,rmse=v['rmse_per10k'],dual=v['dual_probability'],estimable=int(mask.sum()),rank_vs_six=1+int(sum(scores<v['rmse_per10k'])),rmse_mcse=v['rmse_mcse']))
                if m==12:
                    for depth in [20,30,50,100]:
                        result=baseline if depth==k else record(w,ok,depth)
                        areas.append(dict(section=s,method=method,sites=m,k=depth,total_umi=m*depth,footprint_mm2=m*spot_area,**result))
    pd.DataFrame(rows).to_csv(args.out/'Spatial_weighting.csv',index=False)
    pd.DataFrame(areas).to_csv(args.out/'Fixed_area.csv',index=False)
if __name__=='__main__':main()
