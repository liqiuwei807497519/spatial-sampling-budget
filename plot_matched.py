from pathlib import Path
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from design_atlas import policies
from plot_theme import apply_jtm_theme,CMAP_DIV,CMAP_SEQ
apply_jtm_theme()
short={'random':'Random','systematic_grid':'Grid','spatial_only':'Spatial','stratified':'Stratified','stratified_dispersed':'Strata + space','maximin':'Maximin','inhibition':'Inhibition'}
legacy=list(short)[:5];new=list(short)
cols=['#64748B','#C2410C','#0B4F6C','#7E22CE','#0D9488','#4A7BB6','#BB6C28']
def label(fig,x,y,s):
    fig.text(x,y,s,va='top',fontweight='bold',fontsize=10,color='#172938')
def fig4(source, result_root, dest, section=None):
    f=plt.figure(figsize=(170/25.4,205/25.4))
    label(f,.028,.984,'A   Twelve sites at 50 UMIs per site')
    frame=pd.read_csv(source)
    if section is None:section=frame.section.sort_values().iloc[0]
    frame=frame.query("section==@section and spatial_cell>=0 and total-Acan-Col2a1>=50").reset_index(drop=True)
    for i,method in enumerate(legacy):
        ax=f.add_axes([.025+.192*i,.827,.185,.115])
        ids,_,_=policies(frame,12,1,124+i,method)
        colors=np.where(frame.stratum.eq('dark'),'#83B6C2','#DED7C9')
        ax.scatter(frame.x,frame.y,s=1.9,c=colors,lw=0,rasterized=True)
        v=frame.iloc[ids[0]];ax.scatter(v.x,v.y,s=13,c=cols[i],ec='white',lw=.45)
        ax.set(aspect='equal');ax.invert_yaxis();ax.axis('off');ax.set_title(short[method],fontsize=7.2,pad=2)
    legend=[Line2D([],[],marker='o',ls='',color=c,label=l,markersize=4) for c,l in [('#83B6C2','Dark H&E'),('#DED7C9','Light H&E')]]
    f.legend(handles=legend,ncol=2,loc='center',bbox_to_anchor=(.50,.798),handletextpad=.3,columnspacing=1.5)
    raw=pd.read_csv(result_root/'original_five.csv')
    wide=raw.pivot(index=['section','gene'],columns='method',values='rmse')
    gains=100*(1-wide.div(wide.random.replace(0,np.nan),axis=0))
    comp=gains.xs('Acan+Col2a1',level='gene')
    label(f,.028,.765,'B   Paired changes by section')
    ax=f.add_axes([.102,.538,.305,.176])
    for row in comp[legacy].to_numpy():
        ax.plot(range(5),row,color='#A8B9C3',lw=.7,alpha=.85,marker='o',ms=2)
    ax.plot(range(5),comp[legacy].median(),color='#7E22CE',lw=2,marker='D',ms=4,label='Median')
    ax.axhline(0,color='#64748B',lw=.55,ls='--')
    ax.set(xticks=range(5),xticklabels=['Random','Grid','Spatial','Strata','Strata +\nspace'],ylabel='RMSE reduction (%)')
    ax.tick_params(axis='x',rotation=27);ax.legend(loc='lower right')
    label(f,.485,.765,'C   Median gain by RNA target')
    genes=['Acan','Col2a1','Col1a1','Comp','Bgn','Dcn','Fmod','Col9a1','Sox9','Krt19','Mmp3','Mmp13','Adamts5']
    mat=gains.groupby('gene').median().loc[genes,legacy[1:]]
    ax=f.add_axes([.603,.486,.315,.234]);im=ax.imshow(mat,cmap=CMAP_DIV.reversed(),norm=TwoSlopeNorm(vmin=-12,vcenter=0,vmax=12),aspect='auto')
    ax.set(yticks=range(13),yticklabels=genes,xticks=range(4),xticklabels=['Grid','Spatial','Strata','Strata +\nspace'])
    ax.tick_params(length=0,pad=2,labelsize=6.8)
    for tick in ax.get_yticklabels():tick.set_fontstyle('italic')
    cb=f.colorbar(im,cax=f.add_axes([.935,.508,.012,.19]),ticks=[-12,0,12]);cb.ax.tick_params(labelsize=6.5)
    cb.ax.set_title('%',fontsize=6.8,pad=4)
    label(f,.028,.443,'D   Trade-offs at 600 UMIs')
    label(f,.675,.443,'E   Twelve-site footprint')
    pareto=pd.read_csv(result_root/'Pareto_medians.csv')
    a=f.add_axes([.101,.179,.495,.208])
    for method,color in zip(new,cols):
        t=pareto[pareto.method==method]
        a.scatter(t.rmse_per10k,100*t.dual_probability,s=20*t.area_vs_six,color=color,alpha=.84,edgecolors='white',linewidth=.55)
        h=t[t.median_front]
        a.scatter(h.rmse_per10k,100*h.dual_probability,s=20*h.area_vs_six+17,facecolors='none',edgecolors='#1E293B',linewidth=.8)
    a.set(xlabel='Regional contrast RMSE (per 10,000)',ylabel='Request-level dual detection (%)')
    a.grid(axis='both',lw=.38,color='#E5EBEF');a.set_axisbelow(True)
    keys=[Line2D([],[],marker='o',ms=4.5,color=c,ls='',label=short[m]) for m,c in zip(new,cols)]
    f.legend(handles=keys,ncol=4,bbox_to_anchor=(.51,.075),loc='center',fontsize=6.5,handletextpad=.22,columnspacing=.9)
    f.text(.105,.128,'Footprint',fontsize=6.8,va='center')
    for x,m in zip([.235,.341,.446,.555],[6,12,20,30]):
        art=Line2D([x],[.128],transform=f.transFigure,ls='',marker='o',mfc='none',mec='#475569',ms=(20*m/6)**.5,mew=.7)
        f.add_artist(art);f.text(x+.018,.128,f'{m} sites',fontsize=6.5,va='center')
    ar=pd.read_csv(result_root/'Fixed_area.csv').query("method=='stratified_dispersed'").groupby('k').median(numeric_only=True)
    for y,col,lbl,c in [(.293,'rmse_per10k','RMSE','#0B4F6C'),(.157,'dual_probability','Dual detection (%)','#7E22CE')]:
        a=f.add_axes([.747,y,.205,.080])
        v=ar[col]*(100 if col=='dual_probability' else 1)
        a.plot(np.arange(len(v)),v,color=c,marker='o',lw=1.4,ms=3.5)
        a.set(xticks=range(4),xticklabels=[str(i) for i in v.index],ylabel=lbl)
        a.tick_params(labelsize=6.5,pad=1);a.yaxis.label.set_size(6.8)
        if y<.2:a.set_xlabel('UMIs per site',fontsize=6.8,labelpad=1)
    f.text(.045,.022,'Ring: nondominated across error, recovery and sampled spot area. E: stratified dispersion; medians across eight sections.',fontsize=6.5)
    for ext in ('pdf','svg','png'):f.savefig(dest/('Fig4.'+ext),dpi=330 if ext!='png' else 300)
    plt.close(f)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--input',type=Path,default=Path(__file__).parent/'synthetic_pilot.csv');ap.add_argument('--section');args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    fig4(args.input,args.results,args.out,args.section)
