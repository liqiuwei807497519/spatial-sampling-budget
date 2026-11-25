"""Rebuild the supplied rat site frame directly from the eight native Spatial assay objects."""
import argparse,subprocess,tempfile
from pathlib import Path
import numpy as np,pandas as pd
from PIL import Image
from scipy.spatial import cKDTree
from design_atlas import GENES
def main():
    ap=argparse.ArgumentParser();ap.add_argument('rds_directory',type=Path);ap.add_argument('--rscript',default='Rscript');ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='ivd_spots_') as raw:
        tmp=Path(raw)
        encode=lambda x:str(Path(x).resolve()).encode().hex()
        subprocess.run([args.rscript,'--vanilla',str(Path(__file__).with_suffix('.R')),encode(args.rds_directory),encode(tmp)],check=True)
        rows=[]
        order=['Control-1','Control-2','Day_1-1','Day_1-2','Day_3-1','Day_3-2','Day_28-1','Day_28-2']
        for sec in order:
            d=pd.read_csv(tmp/(sec+'.csv')).sort_values('barcode').reset_index(drop=True)
            xy=d[['full_x','full_y']].to_numpy()*d.lowres_scale.iloc[0]
            pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);radius=.275*pitch
            im=np.asarray(Image.open(tmp/(sec+'.png')).convert('RGB'),float)
            od=-np.log((im+1)/256).mean(2);val=np.full(len(d),np.nan)
            for i,(x,y) in enumerate(xy):
                if x-radius<0 or y-radius<0 or x+radius>=im.shape[1] or y+radius>=im.shape[0]:continue
                l,t=np.ceil([x-radius,y-radius]).astype(int);r,b=np.floor([x+radius,y+radius]).astype(int)
                yy,xx=np.mgrid[t:b+1,l:r+1];use=(yy-y)**2+(xx-x)**2<=radius**2
                val[i]=od[yy[use],xx[use]].mean()
            ok=np.flatnonzero(np.isfinite(val));key=d.barcode.to_numpy()
            lo,buf,hi=np.array_split(ok[np.lexsort((key[ok],val[ok]))],3)
            lab=np.full(len(d),'outside',object);cells=np.full(len(d),-1)
            lab[buf]='buffer'
            for j,(group,label) in enumerate([(lo,'light'),(hi,'dark')]):
                lab[group]=label;order_x=group[np.lexsort((key[group],xy[group,0]))]
                for h,half in enumerate(np.array_split(order_x,2)):
                    for s,sites in enumerate(np.array_split(half[np.lexsort((key[half],xy[half,1]))],3)):cells[sites]=6*j+3*h+s
            frame=pd.DataFrame(dict(section=sec,site=key,x=xy[:,0],y=xy[:,1],mean_rgb_od=val,stratum=lab,spatial_cell=cells,total=d.total))
            rows.append(frame.join(d[GENES]))
        tab=pd.concat(rows,ignore_index=True)
        tab.to_csv(args.out,index=False,compression={'method':'gzip','mtime':0})
        print('Exported',len(tab),'registered positions.')
if __name__=='__main__':main()
