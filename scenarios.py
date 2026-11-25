"""Compare measured site-selection designs on a spatial signal with known region contrasts."""
from pathlib import Path
import itertools, argparse
import numpy as np
import pandas as pd
from design_atlas import METHODS,policies
from sampling_design import partition

def main(out,source):
    rng=np.random.default_rng(11857)
    side=32
    y,x=np.mgrid[:side,:side]
    freq=np.fft.fftfreq(side)*2*np.pi
    fy,fx=np.meshgrid(freq,freq,indexing='ij')
    rows=[]
    for shape,fraction in itertools.product(['square','crescent'],[.3,.5]):
        mask=np.ones_like(x,dtype=bool) if shape=='square' else (((x-16)**2+(y-15)**2<225)&((x-23)**2+(y-9)**2>40))
        xy=np.column_stack([x[mask],y[mask]])
        f=pd.DataFrame(xy,columns=['x','y']);f['site']=range(len(f))
        f['stratum']=np.where(f.x<np.quantile(f.x,fraction),'dark','light')
        designs={method:policies(f,12,300,521+j,method) for j,method in enumerate(METHODS)}
        dark=f.stratum.eq('dark').to_numpy()
        for ell,anis,signal,hot in itertools.product([2,8],[1,3],[.25,1.],[False,True]):
            # Periodic lattice Matérn-type spectrum, nu=1.5. Independent fields;
            # each of 300 layouts is paired across scenarios at a fixed frame.
            multiplier=(1+(ell*fx)**2+(ell*anis*fy)**2)**(-1.25)
            white=rng.normal(size=(300,side,side))
            field=np.fft.ifft2(np.fft.fft2(white,axes=(-2,-1))*multiplier,axes=(-2,-1)).real
            field=field[:,mask]
            field=(field-field.mean(1)[:,None])/np.maximum(1e-9,field.std(1)[:,None])
            target=signal*field+(1-signal**2)**.5*rng.normal(size=field.shape)+dark*.4
            if hot:target+=2*np.exp(-((xy[:,0]-10)**2+(xy[:,1]-12)**2)/18)
            truth=target[:,dark].mean(1)-target[:,~dark].mean(1)
            for method,(ids,w,ok) in designs.items():
                estimated=(w*target[np.arange(300)[:,None],ids]).sum(1)
                err=(estimated[ok]-truth[ok])**2
                rows.append(dict(shape=shape,dark_fraction=dark.mean(),length=ell,anisotropy=anis,spatial_coefficient=signal,
                                 hotspot=hot,method=method,estimable=ok.sum(),fields=300,rmse=np.sqrt(err.mean())))
    pd.DataFrame(rows).to_csv(out/'simulation_scenarios.csv',index=False)
    # How many geographical cells per image stratum at the same six-site quota?
    data=pd.read_csv(source)
    outrows=[]
    from design_atlas import evaluate
    for i,(sec,frame) in enumerate(data.query('spatial_cell>=0').groupby('section')):
        f=frame[frame.total-frame.Acan-frame.Col2a1>=50].reset_index(drop=True)
        xy=f[['x','y']].to_numpy();population=[np.flatnonzero(f.stratum.to_numpy()==s) for s in ('dark','light')]
        for k in (1,2,3,6):
            bins=[b for g in population for b in partition(g,xy,k)]
            per=6//k
            gen=np.random.default_rng(32100+i*21+k)
            ids=np.column_stack([np.stack([gen.choice(b,per,replace=False) for _ in range(1500)]) for b in bins])
            mass=np.repeat([len(b)/per for b in bins],per)
            w=np.tile(np.r_[mass[:6]/mass[:6].sum(),-mass[6:]/mass[6:].sum()],(len(ids),1))
            val=evaluate(f,ids,w,np.ones(len(w),bool),50,{'primary':('Acan','Col2a1')})[0]
            outrows.append(dict(section=sec,cells_per_stratum=k,**val))
    pd.DataFrame(outrows).to_csv(out/'partition_sensitivity.csv',index=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument("--input",type=Path,default=Path(__file__).parent/"synthetic_pilot.csv");args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True);main(args.out,args.input)
