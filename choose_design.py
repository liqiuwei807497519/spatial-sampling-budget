"""Compare predeclared layouts on a registered pilot, under molecular and footprint constraints."""
import argparse,math
from pathlib import Path
import pandas as pd
from design_atlas import METHODS,policies,evaluate
from resource_analysis import nondominated
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True);ap.add_argument('--section',required=True)
    ap.add_argument('--targets',nargs=2,required=True);ap.add_argument('--budget',type=int,required=True)
    ap.add_argument('--sites',type=int,nargs='+',required=True);ap.add_argument('--max-sites',type=int)
    ap.add_argument('--spot-diameter-um',required=True,type=float);ap.add_argument('--draws',type=int,default=2000)
    ap.add_argument('--seed',type=int,default=42);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();a,b=args.targets
    if a==b or args.spot_diameter_um<=0 or args.draws<2 or args.budget<=0:ap.error('Two distinct targets and positive budgets are required.')
    m_values=sorted(m for m in set(args.sites) if args.max_sites is None or m<=args.max_sites)
    if not m_values or any(m<2 or m%2 or args.budget%m for m in m_values):ap.error('Choose positive even site counts dividing the budget.')
    frame=pd.read_csv(args.input);frame=frame[(frame.section==args.section)&frame.stratum.isin(['dark','light'])].copy()
    count=frame[[a,b,'total']]
    if count.isna().any().any() or (count<0).any().any() or ((count%1)!=0).any().any() or (count[a]+count[b]>count.total).any():ap.error('Invalid target or total counts.')
    frame=frame[frame.total-frame[a]-frame[b]>=args.budget//min(m_values)].reset_index(drop=True)
    if set(frame.stratum)!= {'dark','light'} or frame.stratum.value_counts().min()<max(m_values)//2:ap.error('Insufficient sites in the two pilot strata.')
    rows=[]
    for m in m_values:
        k=args.budget//m
        for j,method in enumerate(METHODS):
            sel,w,ok=policies(frame,m,args.draws,args.seed+13*m+j,method)
            v=evaluate(frame,sel,w,ok,k,{'specified':(a,b)})[0]
            rows.append(dict(method=method,sites=m,k=k,total_umi=args.budget,footprint_mm2=m*math.pi*(args.spot_diameter_um/2000)**2,**v))
    z=pd.DataFrame(rows)
    values=z[['rmse_per10k','dual_probability','footprint_mm2']].to_numpy(copy=True);values[:,1]*=-1
    z['on_pilot_front']=nondominated(values)
    z['lowest_pilot_error']=z.rmse_per10k==z.rmse_per10k.min()
    z['highest_pilot_detection']=z.dual_probability==z.dual_probability.max()
    args.out.parent.mkdir(parents=True,exist_ok=True);z.to_csv(args.out,index=False)
    print(z[z.on_pilot_front][['method','sites','k','footprint_mm2','rmse_per10k','dual_probability']].to_string(index=False))
if __name__=='__main__':main()
