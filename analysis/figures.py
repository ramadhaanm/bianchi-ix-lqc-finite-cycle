"""Figures and generated LaTeX tables from completed v3 reanalysis only."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT/'paper/figs'
GEN = ROOT/'paper/generated'
COLORS = {'K':'#2368ac','A':'#b74c26'}


def save(fig,name):
    FIG.mkdir(exist_ok=True)
    for ext in ('pdf','png'):
        fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',dpi=200)
    plt.close(fig)


def make_figures(summary,runs,basin,include_long=True):
    plt.rcParams.update({'font.size':8,'pdf.fonttype':42,'ps.fonttype':42,
        'axes.spines.top':False,'axes.spines.right':False,
        'axes.labelsize':8,'axes.titlesize':8,'legend.fontsize':7})
    fig,axs=plt.subplots(1,4,figsize=(7.1,2.0),sharey=True)
    classical=np.loadtxt(ROOT/'results/cl_basins_fig2_200_d64.csv',delimiter=',',comments='#',dtype=int)
    cmap=ListedColormap(['#9ec5f4','#d85045','black','#999999','white'])
    axs[0].imshow(classical,origin='lower',extent=[1,2,1,2],cmap=cmap,vmin=-1.5,vmax=3.5,interpolation='nearest')
    axs[0].set_title('Classical GR')
    for ax,e,title in zip(axs[1:],[1,3,6],['K: first bounce','K: first recollapse','K: second bounce']):
        grid=np.zeros((basin.N,basin.N),int)
        for k in basin.ids:
            r=basin.rows[k,'base']
            a=basin.event(k,'base',e)
            lab=int(np.argmin(a[:6:2]+a[1:6:2]))+1 if basin.valid(k,'base',e) else (-1 if r['reason']=='physical' else 0)
            grid.flat[k-1]=lab
        ax.imshow(grid,origin='lower',extent=[1,2,1,2],cmap=cmap,vmin=-1.5,vmax=3.5,interpolation='nearest')
        ax.set_title(title)
    for ax in axs:
        ax.set_xlabel('$u$'); ax.set_xticks([1.2,1.5,1.8]); ax.set_yticks([1.2,1.5,1.8])
    axs[0].set_ylabel('$v$')
    fig.text(.5,-.04,'Black / grey / white: axis 1 / 2 / 3; blue: physical exit; red: numerical or unresolved.',ha='center',fontsize=6.5)
    save(fig,'fig1_basins')

    fig,axs=plt.subplots(1,2,figsize=(7.1,2.5),sharey=True)
    palette=['#ce5b32','#21866b','#4e63b0']
    for ax,q in zip(axs,('K','A')):
        for c,color in zip(range(1,4),palette):
            rows=summary['labels'][q]['events'][f'B{c}']['series'][1:]
            rows=[r for r in rows if r['n']]
            positive=[r for r in rows if r['flips']>0]
            if positive:
                x=np.array([r['eps'] for r in positive]); y=np.array([r['f'] for r in positive])
                lo=np.array([r['ci'][0] for r in positive]); hi=np.array([r['ci'][1] for r in positive])
                ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt='o-',ms=3,lw=.8,color=color,label=f'Bounce {c}')
            zero=[r for r in rows if r['flips']==0]
            if zero:
                ax.scatter([r['eps'] for r in zero],[r['ci'][1] for r in zero],marker='v',s=16,facecolors='none',edgecolors=color)
        ax.axhline(.1,lw=.6,ls='--',color='.5')
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-3,1)
        ax.set_xlabel(r'Initial separation $\varepsilon$'); ax.set_title(rf'{q} quantization, $\Lambda=20$')
        ax.grid(alpha=.2)
    axs[0].set_ylabel('Label disagreement fraction')
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.15),ncol=3)
    save(fig,'fig2_uncertainty')

    order=['turnaround','contraction','bounce','expansion']
    labels=['Turnaround\n'+r'($\Delta\ln V=5$)','Deep\ncontraction','Bounce\nwindow','Expansion']
    fig,axs=plt.subplots(1,2,figsize=(7.1,2.6),sharey=True)
    for ax,q in zip(axs,('K','A')):
        for j,L in enumerate((10,20,30)):
            d=summary['phases'][f'{q}{L}']['cutoffs']['1']
            vals=[d[p]['median'] for p in order]
            err=np.array([[d[p]['median']-d[p]['lo'] for p in order],[d[p]['hi']-d[p]['median'] for p in order]])
            x=np.arange(4)+(j-1)*.27
            ax.bar(x,vals,width=.25,color=['#a9c9e8','#508dc3','#1c5186'][j],label=rf'$\Lambda={L}$')
            ax.errorbar(x,vals,yerr=err,fmt='none',color='.2',lw=.8,capsize=2)
        ax.set_xticks(range(4),labels); ax.tick_params(axis='x',labelsize=7)
        ax.set_title(f'{q} quantization'); ax.grid(axis='y',alpha=.2)
    axs[0].set_ylabel('Median log$_{10}$ amplification\n(95% bootstrap confidence interval)')
    axs[1].legend()
    save(fig,'fig3_phases')

    fig,ax=plt.subplots(figsize=(3.4,2.5))
    for q,m in [('K','o'),('A','s')]:
        d=[summary['phases'][f'{q}{L}'] for L in (10,20,30)]
        x=[r['cohort_volume_depth']['median'] for r in d]
        xerr=[[r['cohort_volume_depth']['median']-r['cohort_volume_depth']['lo'] for r in d],
              [r['cohort_volume_depth']['hi']-r['cohort_volume_depth']['median'] for r in d]]
        y=[r['cutoffs']['1']['total']['median'] for r in d]
        err=[[r['cutoffs']['1']['total']['median']-r['cutoffs']['1']['total']['lo'] for r in d],
             [r['cutoffs']['1']['total']['hi']-r['cutoffs']['1']['total']['median'] for r in d]]
        ax.errorbar(x,y,xerr=xerr,yerr=err,marker=m,color=COLORS[q],label=q,capsize=3)
    ax.set_xlabel('Measured initial volume contraction $N_V$')
    ax.set_ylabel('Median cycle log$_{10}$ amplification')
    ax.grid(alpha=.2); ax.legend()
    sec=ax.secondary_xaxis('top',functions=(lambda x:x/3,lambda x:3*x))
    sec.set_xlabel('Mean-scale-factor contraction $N_a=N_V/3$')
    save(fig,'fig4_rate_vs_depth')

    if not include_long:
        return

    fig,axs=plt.subplots(1,2,figsize=(7.1,2.5),sharey=True)
    for ax,ev in zip(axs,('B','R')):
        for q,dx in [('K',-.06),('A',.06)]:
            st=summary['growth'][q]['steps'][ev]
            x=np.arange(1,6)+dx
            vals={key:np.array([s[key] if s[key] is not None else np.nan for s in st])
                  for key in ('median','lo','hi')}
            y=vals['median']
            ax.errorbar(x,y,yerr=[y-vals['lo'],vals['hi']-y],marker='o' if q=='K' else 's',label=q,color=COLORS[q],capsize=3)
            for xx,yy,s in zip(x,y,st):
                if np.isfinite(yy):
                    ax.annotate(str(s['n']),(xx,yy),xytext=(-7 if q=='K' else 7,7 if q=='K' else -12),
                                textcoords='offset points',fontsize=6,ha='center',color=COLORS[q],
                                bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.2))
        ax.set_xticks(range(1,6),[f'{c}\u2013{c+1}' for c in range(1,6)])
        ax.set_xlabel('Consecutive cycles'); ax.set_title('Bounces' if ev=='B' else 'Recollapses'); ax.grid(alpha=.2)
    axs[0].set_ylabel('Median log$_{10}$ amplification per cycle')
    axs[1].legend(); save(fig,'fig5_longrun')

def make_tables(summary,write_report=True):
    from reanalyze import fmt
    GEN.mkdir(exist_ok=True)
    phase=[]; cut=[]; survival=[]; depths=[]; threshold=[]; exits=[]; labelrows=[]; quality=[]; long=[]; longsurvival=[]
    for q in ('K','A'):
        for L in (10,20,30):
            d=summary['phases'][f'{q}{L}']; z=d['cutoffs']['1']
            phase.append(f"{q} & {L} & {z['resolved']}/{z['eligible']} & {fmt(z['bounce'])} & {fmt(z['total'])} " + r'\\')
            survival.append(f"{q} & {L} & "+' & '.join(str(d['survival'][str(c)]) for c in (1,2,3))+r' \\')
            depths.append(f"{q} & {L} & {fmt(d['volume_depth'],False)} & {fmt(d['scale_depth'],False)} "+r'\\')
            exits.append(f"{q} & {L} & "+' & '.join(str(d['stop'][s]) for s in ('complete','physical','constraint','budget','integration'))+r' \\')
        for cutoff in (1,3,10):
            z=summary['phases'][f'{q}20']['cutoffs'][str(cutoff)]
            cut.append(f"{q} & {cutoff} & {z['resolved']}/{z['eligible']} & {fmt(z['bounce'])} & {fmt(z['total'])} "+r'\\')
        for name,z in summary['phases'][q+'20']['quality_sensitivity'].items():
            quality.append(f"{q} & {name} & {z['n']} & {fmt(z['bounce'])} & {fmt(z['total'])} "+r'\\')
        g=summary['growth'][q]
        for ev in ('B','R'):
            for cutoff in (1,3,10):
                z=g['cutoffs'][str(cutoff)][ev]
                long.append(f"{q} & {ev} & {cutoff} & {z['n']}/{z['eligible']} & {fmt(z)} "+r'\\')
        longsurvival.append(q+' & '+' & '.join(str(g['survival'][str(c)]) for c in range(1,7))+' & '+
            ' & '.join(str(g['stop'][s]) for s in ('complete','physical','constraint','budget','integration'))+r' \\')
        for onset in (.3,.5,.7):
            z=summary['matched_thresholds'][q][str(onset)]
            threshold.append(f"{q} & {onset} & {z['resolved']} & {fmt(z['bounce'])} & {fmt(z['total'])} "+r'\\')
        for ev in ('B1','B2','B3','R1','R2','R3'):
            z=summary['labels'][q]['events'][ev]
            def sci(x):
                exp=int(np.floor(np.log10(x))); return rf'{x/10**exp:.0f}\times10^{{{exp}}}'
            bracket='; '.join('$['+sci(a)+','+sci(b)+']$' for a,b in z['brackets']) or 'unbracketed'
            labelrows.append(f"{q} & ${ev[0]}_{ev[1]}$ & {bracket} & {z['series'][0]['flips']}/{z['series'][0]['n']} "+r'\\')
    for name,lines in [('phase_rows',phase),('cutoff_rows',cut),('survival_rows',survival),('depth_rows',depths),('threshold_rows',threshold),('exit_rows',exits),('label_rows',labelrows),('quality_rows',quality),('long_rows',long),('long_survival_rows',longsurvival)]:
        (GEN/f'{name}.tex').write_text('\n'.join(lines)+'\n')
    # Headline macros are derived, not manually transcribed.
    macro=[]
    for q in ('K','A'):
        d=summary['phases'][f'{q}20']['cutoffs']['1']
        for name,key in [('Bounce','bounce'),('Cycle','total')]:
            macro.append('\\newcommand{\\'+q+name+'}{'+fmt(d[key],False)+'}')
        for ev in ('B','R'):
            z=summary['growth'][q]['steps'][ev+'_cohort']
            macro.append('\\newcommand{\\'+q+'Long'+ev+'}{'+fmt(z,False)+'}')
    (GEN/'numbers.tex').write_text('\n'.join(macro)+'\n')
    report=['Corrected reanalysis: all quoted intervals are 95% percentile-bootstrap intervals of the median.','',
        'Pathwise threshold and constraint exits are irreversible. Event survival refers to reaching R_n before an exit.',
        'Counts of complete/physical/constraint/budget/integration below refer to the whole requested run; phase runs also attempt P_3.','']
    for name,d in summary['phases'].items():
        z=d['cutoffs']['1']
        report += [f"{name}: stop reasons {d['stop']}; R1/R2/R3 survival {d['survival']}; resolved/eligible {z['resolved']}/{z['eligible']}; bounce {fmt(z['bounce'])}; cycle {fmt(z['total'])}."]
    report += ['','Quality-selection sensitivity (supplementary; primary thresholds determine the manuscript statistics):']
    for q in ('K','A'):
        for name,z in summary['phases'][q+'20']['quality_sensitivity'].items():
            report += [f"{q}, {name}: n={z['n']}, bounce {fmt(z['bounce'])}, cycle {fmt(z['total'])}."]
    report += ['','Six-cycle rates use a fixed cohort passing quality controls at every respective event.']
    for q,d in summary['growth'].items():
        report += [f"{q}: stop {d['stop']}; survival {d['survival']}; B {fmt(d['steps']['B_cohort'])}, n={d['steps']['B_cohort']['n']}; R {fmt(d['steps']['R_cohort'])}, n={d['steps']['R_cohort']['n']}."]
        for cutoff in (3,10):
            z=d['cutoffs'][str(cutoff)]
            report += [f"{q}, six-cycle cutoff {cutoff}: B {fmt(z['B'])}, n={z['B']['n']}; R {fmt(z['R'])}, n={z['R']['n']}."]
    if write_report:
        (ROOT/'results/revised/RESULTS.md').write_text('\n\n'.join(report)+'\n')
