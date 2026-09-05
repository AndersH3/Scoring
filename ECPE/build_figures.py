"""Rebuild vector figures and exact-data LaTeX tables for the ECPE report."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
D, S, F, T = [ROOT / p for p in ('data', 'supplementary', 'figures', 'tables')]
for p in (F,T): p.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9, 'axes.spines.top':False,
                     'axes.spines.right':False, 'pdf.fonttype':42, 'axes.labelsize':9,
                     'legend.fontsize':8, 'axes.titlesize':10})
blue, orange, grey = '#246482', '#B9512D', '#77838C'
def read(name, folder=D): return pd.read_csv(folder / (name+'.csv'))
def save(fig,name):
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)
def table(name, headers, rows, spec=None):
    n=len(headers); spec=spec or ('l'+'r'*(n-1))
    body='\\begin{longtblr}[caption={'+name.replace('_',' ').capitalize()+'}]{colspec={'+spec+'},rowhead=1,row{1}={font=\\bfseries,bg=lightblue},hlines,rowsep=2pt}\n'
    body+=' & '.join(headers)+' \\\\\n'
    body+='\n'.join(' & '.join(map(str,r))+' \\\\' for r in rows)
    body+='\n\\end{longtblr}\n'
    (T/(name+'.tex')).write_text(body)

w=read('item_weights_full'); p=read('person_scores'); b=read('bootstrap_item_intervals'); q=read('ecpe_qmatrix')
fig,axs=plt.subplots(1,2,figsize=(7.2,3.05),layout='constrained')
axs[0].scatter(p.number_correct,p.full_weighted,s=5,alpha=.13,c=blue,rasterized=True)
g=read('score_by_number_correct');axs[0].plot(g.number_correct,g['mean'],c=orange,label='Mean within raw total')
axs[0].plot([0,28],[0,1],'--',c=grey,label='Raw proportion')
axs[0].set(xlabel='Number correct',ylabel='Weighted score',xlim=(0,28),ylim=(0,1));axs[0].legend(loc='upper left')
sc=axs[1].scatter(w.proportion_correct,100*w.share,c=w.Dxy,cmap='viridis',s=40)
fig.colorbar(sc,ax=axs[1],label='Item-rest D',shrink=.8)
for _,r in w.nlargest(3,'share').iterrows():axs[1].annotate(r['item'],(r.proportion_correct,100*r.share),xytext=(4,3),textcoords='offset points',fontsize=8)
axs[1].axhline(100/28,ls='--',c=grey)
axs[1].set(xlabel='Proportion correct',ylabel='Item share (%)',xlim=(.39,.95))
save(fig,'baseline')

bs=b[b.quantity=='share'].set_index('item').loc[w.item]
fig,ax=plt.subplots(figsize=(7.2,3.5),layout='constrained')
xx=np.arange(28);ax.errorbar(xx,100*bs.estimate,yerr=np.stack((100*(bs.estimate-bs.p025),100*(bs.p975-bs.estimate))),fmt='o',ms=3,c=blue,capsize=2)
ax.axhline(100/28,ls='--',c=grey);ax.set_xticks(xx,w.item,rotation=90)
ax.set(ylabel='Calibration item share (%)',xlabel='Item');save(fig,'bootstrap_weights')

ss=read('sample_size_replicates'); fig,axs=plt.subplots(1,2,figsize=(7.2,3),layout='constrained')
ns=sorted(ss.calibration_n.unique())
for ax,col,scale,label in zip(axs,['rmse','p95_rank_shift_pp'],[100,1],['Audit RMSE (score percentage points)','95th percentile rank shift (pp)']):
    ax.boxplot([scale*ss.loc[ss.calibration_n==n,col] for n in ns],tick_labels=list(map(str,ns)),showfliers=False,patch_artist=True,boxprops={'facecolor':'#D8E7ED'},medianprops={'color':orange})
    ax.set(xlabel='Calibration persons (without replacement)',ylabel=label)
save(fig,'sample_size')

co=read('contamination_replicates'); fig,ax=plt.subplots(figsize=(7.2,3.1),layout='constrained')
for kind,label,color in [('cell_flip','Cell flips',blue),('random_person_50pct','Random response rows',orange),('all_zero_person','All-zero rows','#7761A6'),('all_one_person','All-one rows','#568458')]:
    z=co[co.kind==kind].groupby('rate').rmse
    med=z.median(); low=z.quantile(.025);high=z.quantile(.975)
    ax.plot(100*med.index,100*med.values,'o-',label=label,c=color)
    ax.fill_between(100*med.index,100*low.values,100*high.values,alpha=.12,color=color)
ax.set(xlabel='Requested calibration corruption (%)',ylabel='Audit RMSE (score percentage points)',xticks=[1,5,10]);ax.legend(ncol=2)
save(fig,'corruption')

nm=read('permutation_null_metrics'); sh=read('split_half_consistency')
fig,axs=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
axs[0].boxplot([nm.mean_Dxy_rest,nm.mean_Dxy_total],tick_labels=['Item-rest','Item-total'],patch_artist=True,boxprops={'facecolor':'#D8E7ED'},medianprops={'color':orange});axs[0].axhline(0,c=grey,ls='--');axs[0].set(ylabel='Mean item Somers D',title='Independent-column null')
axs[1].hist(sh.corrected_difference,bins=16,color=blue,alpha=.85);axs[1].axvline(0,c=orange,ls='--');axs[1].set(xlabel='Weighted minus raw correction',ylabel='Random halves',title='Paired split-half differences')
save(fig,'null_reliability')

au=read('augmented_item_stress'); labels=['All zero','All one','Random 1%','Random 50%','Reversed','One duplicate','Five duplicates']
fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),layout='constrained')
for ax,col,mul,xlabel in zip(axs,['extra_share','rmse'],[100,100],['Added items: total score share (%)','Audit RMSE (score percentage points)']):
    ax.barh(labels,mul*au[col],color=blue);ax.invert_yaxis();ax.set(xlabel=xlabel)
save(fig,'augmentation')

nu=read('null_item_replicates',S)
fig,axs=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
rates=sorted(nu.probability.unique())
axs[0].boxplot([100*nu.loc[nu.probability==v,'share'] for v in rates],tick_labels=['0.1%','1%','5%','50%'],showfliers=False,patch_artist=True,boxprops={'facecolor':'#D8E7ED'},medianprops={'color':orange});axs[0].set(xlabel='Independent success probability',ylabel='Added item share (%)')
for prob,color in [(.001,orange),(.01,blue),(.05,'#568458')]:
    vals=np.sort(100*nu.loc[nu.probability==prob,'share'].to_numpy());axs[1].plot(vals,np.arange(1,len(vals)+1)/len(vals),label=f'{100*prob:g}%',c=color)
axs[1].set(xlabel='Added item share (%)',ylabel='Empirical cumulative fraction',ylim=(0,1));axs[1].legend(title='Null probability')
save(fig,'null_added')

rp=read('regression_audit_predictions',S); coeff=read('regression_coefficients',S).iloc[0];xx=np.linspace(0,1,201)
fig,axs=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
axs[0].scatter(rp.q,rp.frozen_weighted_score,s=7,alpha=.2,c=blue,rasterized=True);axs[0].plot(xx,xx+coeff.kappa*xx*(1-xx),c=orange,label='Constrained quadratic');axs[0].plot(xx,xx,'--',c=grey,label='Identity');axs[0].set(xlabel='Raw proportion',ylabel='Frozen weighted score',xlim=(0,1),ylim=(0,1));axs[0].legend()
axs[1].scatter(rp.q,rp.frozen_weighted_score-rp.constrained_prediction,s=8,alpha=.3,c=blue,rasterized=True);axs[1].axhline(0,ls='--',c=grey);axs[1].set(xlabel='Raw proportion',ylabel='Weighted score minus fitted curve')
save(fig,'regression')

fig,axs=plt.subplots(1,2,figsize=(7.2,3.0),layout='constrained');n=2045
pp=np.geomspace(.0001,1,500)
for a in [0,.3,.6,.9]:axs[0].plot(pp,(1-np.log((n*pp+1)/(n+1)))*(1-np.log((n*(1-a)+1)/(n+1))),label=f'A = {a:g}')
axs[0].set(xscale='log',xlabel='Proportion correct',ylabel='Raw item weight');axs[0].legend(ncol=2)
aa=np.linspace(0,1,500)
for nn in [50,500,2045]:axs[1].plot(aa,1-np.log((nn*(1-aa)+1)/(nn+1)),label=f'n = {nn}')
axs[1].set(xlabel='Positive Somers component A',ylabel='Discrimination factor');axs[1].legend()
save(fig,'formula')

rows=[]
for _,r in w.iterrows():
    skills=''.join(str(int(v)) for v in q.set_index('item').loc[r['item']])
    rows.append([r['item'],skills,int(r.correct),f'{r.proportion_correct:.4f}',f'{r.Dxy:.4f}',f'{r.rarity_factor:.4f}',f'{r.discrimination_factor:.4f}',f'{r.raw_weight:.4f}',f'{100*r.share:.3f}'])
table('complete_item_weights',['Item','Q','Correct','$p$','$D$','$F_1$','$F_2$','$w$','Share \\%'],rows)
table('bootstrap_share_intervals',['Item','Estimate \\%','2.5th \\%','97.5th \\%'],[[k,f'{100*r.estimate:.3f}',f'{100*r.p025:.3f}',f'{100*r.p975:.3f}']for k,r in bs.iterrows()])
table('scores_by_raw_total',['Correct','Persons','Minimum','Mean','Maximum','SD'],[[int(r.number_correct),int(r.persons)]+[f'{r[k]:.4f}' if pd.notna(r[k]) else '---' for k in ['min','mean','max','sd']] for _,r in g.iterrows()])
out={'figures':len(list(F.glob('*.pdf'))),'largest_within_total_span':float((g['max']-g['min']).max())}
for gap in range(27,0,-1):
    found=False
    for total in range(5,29-gap):
        a=p[p.number_correct==total];bb=p[p.number_correct==total+gap]
        if len(a) and len(bb) and a.full_weighted.max()>bb.full_weighted.min()+1e-12:
            lo=a.loc[a.full_weighted.idxmax()];hi=bb.loc[bb.full_weighted.idxmin()]
            out['maximum_raw_gap_reversal']={'gap':gap,'lower_raw_person':int(lo.person),'lower_raw':int(lo.number_correct),'lower_raw_weighted':float(lo.full_weighted),'higher_raw_person':int(hi.person),'higher_raw':int(hi.number_correct),'higher_raw_weighted':float(hi.full_weighted)}
            found=True;break
    if found:break
(S/'report_derived_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
