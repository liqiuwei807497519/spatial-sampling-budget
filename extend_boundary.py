from pathlib import Path
import pandas as pd,numpy as np
from scipy.spatial import cKDTree
from design_atlas import evaluate,policies

def main(source,out):
    tab=pd.read_csv(source);rows=[]
    for i,(section,z) in enumerate(tab.groupby('section')):
        z=z[z.stratum.isin(['dark','light','buffer'])].copy()
        # Match the common budget frame and retain the original clear endpoint labels.
        z=z[z.total-z[['Acan','Col2a1','Sox9','Mmp3','Adamts5','Col9a1']].sum(axis=1)>=100].reset_index(drop=True)
        core=z[z.stratum!='buffer'];buf=z[z.stratum=='buffer'].copy()
        xy=core[['x','y']].to_numpy()
        pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1])
        d=[] 
        for s in ['dark','light']:
            d.append(cKDTree(core.loc[core.stratum.eq(s),['x','y']]).query(buf[['x','y']])[0]/pitch)
        distance=np.minimum(d[0],d[1]);distinct=~np.isclose(d[0],d[1],atol=1e-8)
        buf['stratum']=np.where(d[0]<d[1],'dark','light')
        for radius in [1,2]:
            extension=buf[(distance<=radius+1e-7)&distinct]
            f=pd.concat([core,extension],ignore_index=True).drop(columns='spatial_cell')
            for method in ['random','stratified','stratified_dispersed']:
                sel,w,ok=policies(f,12,2000,17890+i*29+radius,method)
                val=evaluate(f,sel,w,ok,50,{'primary':('Acan','Col2a1')})[0]
                rows.append(dict(section=section,outward_pitches=radius,added_sites=len(extension),method=method,**val))
    pd.DataFrame(rows).to_csv(Path(out)/'outward_sensitivity.csv',index=False)
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--input',type=Path,default=Path(__file__).parent/'synthetic_pilot.csv');args=ap.parse_args()
    main(args.input,args.out)
