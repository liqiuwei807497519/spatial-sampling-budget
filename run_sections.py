"""Resume the eight section-level jobs without changing the full-analysis seeds or frame."""
import argparse,sys,subprocess
from pathlib import Path
import pandas as pd
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--input',type=Path,default=Path(__file__).parent/'synthetic_pilot.csv');a=ap.parse_args()
    here=Path(__file__).parent
    sections=sorted(pd.read_csv(a.input,usecols=['section']).section.unique())
    files=['original_five','fixed600','margin_sensitivity','design_lookup']
    for s in sections:
        job=a.out/'sections'/s
        if not (job/'complete').is_file():
            subprocess.run([sys.executable,'-B',str(here/'design_atlas.py'),'--input',str(a.input.resolve()),'--section',s,'--out',str(job)],check=True)
            (job/'complete').write_text('OK\n')
    a.out.mkdir(parents=True,exist_ok=True)
    for fname in files:
        table=pd.concat([pd.read_csv(a.out/'sections'/s/(fname+'.csv')) for s in sections],ignore_index=True)
        table.to_csv(a.out/(fname+'.csv'),index=False)
    f=pd.read_csv(a.out/'fixed600.csv');keys=['section','sites','pair']
    f=f.merge(f[f.method=='random'][keys+['rmse_per10k']].rename(columns={'rmse_per10k':'base'}),on=keys)
    f['gain']=100*(1-f.rmse_per10k/f.base)
    f.groupby(['pair','sites','method'])[['gain','rmse_per10k','dual_probability']].median().to_csv(a.out/'fixed600_medians.csv')
    pd.read_csv(a.out/'sections'/sections[0]/'abundance.csv').to_csv(a.out/'abundance.csv',index=False)
if __name__=='__main__':main()
