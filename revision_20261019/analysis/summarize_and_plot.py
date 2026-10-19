"""Read saved inference results; no fitting, simulation, or selection changes.

Post hoc diagnostics requested by independent review are labelled explicitly.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
MODELS=('M0','M1','M2','M3')
DESIGNS=('A_step_step','B_step_combination')


def write_csv(path,rows):
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,columns);w.writeheader();w.writerows(rows)


def finish(fig,name):
    for extension in ('pdf','svg','png'):
        fig.savefig(HERE/'figures'/f'{name}.{extension}',dpi=300,bbox_inches='tight')
    plt.close(fig)


def main():
    fits=sum([json.loads((HERE/m/'fits.json').read_text('utf-8')) for m in MODELS],[])
    summary=json.loads((HERE/'summary.json').read_text('utf-8'))
    selected=[r for r in fits if r['selected']]
    failure_rows=[]
    for r in fits:
        if not r['success']:
            failure_rows.append({k:r[k] for k in ('generator','case','sigma','replicate','design','candidate')} |
                                {'true_family':r['candidate']==r['generator'],
                                 'statuses':';'.join(str(s.get('status','error')) for s in r['starts']),
                                 'messages':' | '.join(s['message'] for s in r['starts'])})
    summary['candidate_failure_details']=failure_rows
    summary['true_family_failed_fit_count']=sum(r['true_family'] for r in failure_rows)
    lookup={(r['generator'],r['case'],r['sigma'],r['replicate'],r['design'],r['candidate']):r for r in fits}
    paired=[]
    for scope in ('all_assigned','true_family_converged_both','all_candidates_converged_both'):
        for sigma in (.01,.03):
            keys=sorted(set((r['generator'],r['case'],r['sigma'],r['replicate']) for r in selected if r['sigma']==sigma))
            if scope=='true_family_converged_both':
                keys=[k for k in keys if all(lookup[(*k,d,k[0])]['success'] for d in DESIGNS)]
            elif scope=='all_candidates_converged_both':
                keys=[k for k in keys if all(lookup[(*k,d,m)]['success'] for d in DESIGNS for m in MODELS)]
            choices={d:[next(r for m in MODELS if (r:=lookup[(*k,d,m)])['selected']) for k in keys] for d in DESIGNS}
            a,b=choices[DESIGNS[0]],choices[DESIGNS[1]]
            row=dict(scope=scope,post_hoc=scope!='all_assigned',sigma=sigma,paired_datasets=len(keys),
                     correct_A=sum(r['candidate']==r['generator'] for r in a),correct_B=sum(r['candidate']==r['generator'] for r in b),
                     improved_family_selection=sum(x['candidate']!=x['generator'] and y['candidate']==y['generator'] for x,y in zip(a,b)),
                     worsened_family_selection=sum(x['candidate']==x['generator'] and y['candidate']!=y['generator'] for x,y in zip(a,b)))
            for schedule in ('altered_pulse','full_withdrawal'):
                av=np.array([r[schedule+'_clean_rmse'] for r in a]);bv=np.array([r[schedule+'_clean_rmse'] for r in b])
                row.update({schedule+'_clean_rmse_A_median':float(np.median(av)),schedule+'_clean_rmse_B_median':float(np.median(bv)),
                            schedule+'_paired_B_minus_A_median':float(np.median(bv-av)),schedule+'_B_lower_count':int(sum(bv<av))})
            paired.append(row)
    summary['paired_sensitivity_post_hoc']=paired
    m1_paired=[]
    for sigma in (.01,.03):
        keys=sorted(set((r['generator'],r['case'],r['sigma'],r['replicate']) for r in selected if r['sigma']==sigma and r['generator']=='M1'))
        keys=[k for k in keys if all(lookup[(*k,d,'M1')]['success'] for d in DESIGNS)]
        for schedule in ('altered_pulse','full_withdrawal'):
            a=[lookup[(*k,DESIGNS[0],'M1')][schedule+'_clean_rmse'] for k in keys]
            b=[lookup[(*k,DESIGNS[1],'M1')][schedule+'_clean_rmse'] for k in keys]
            m1_paired.append(dict(post_hoc=True,generator='M1',candidate_fixed='M1',sigma=sigma,schedule=schedule,paired_successful_fits=len(keys),
                                  median_A=float(np.median(a)),median_B=float(np.median(b)),B_lower_count=int(sum(y<x for x,y in zip(a,b)))))
    summary['post_hoc_M1_true_family_paired_prediction']=m1_paired
    summary['interpretation_limits']=[
        'Observed design gains combine statistical discrimination and optimizer reliability; five true-family failures occurred under A, none under B.',
        'Parameter-recovery errors are conditional on successful true-family fitting; denominators must be reported.',
        'The primary small synthetic design is not an estimate of clinical mechanism recovery or a population acceptance probability.',
        'Frozen-M4 a and b and observation map are fixed known; this is not evolving-M4 parameter inference.',
        'A singleton .025-spaced profile-grid set is grid-unresolved, not zero uncertainty. The surrounding bracket can be .05 wide.',
        'The step-only upper profile endpoint near .775 is imposed by r<=4 through r(1-d), not information identifying d.',
        'Historic raw field diagnostic_noise_free_truth_rmse=0 is a generating-map identity only; no noiseless parameter-recovery experiment was performed. It is excluded from all result tables and figures.',
    ]
    for row in summary['parameter_recovery']:
        row.update(assigned_n=12,failed_n=12-row['n'],conditioning='successful_true_family_fits_only')
    prof=json.loads((HERE/'frozen_M4/fits.json').read_text('utf-8'))
    profgrid=list(csv.DictReader((HERE/'frozen_M4/profiles.csv').open(encoding='utf-8')))
    flagged=[]
    for row in prof:
        match=[r for r in profgrid if r['design']==row['design'] and float(r['sigma'])==row['sigma'] and int(r['replicate'])==row['replicate']]
        accepted=[float(r['d']) for r in match if r['in_profile_set']=='True']
        singleton=len(accepted)==1
        flagged.append({k:row[k] for k in ('design','sigma','replicate','r','d','profile_set_min','profile_set_max','profile_set_width')} |
                       dict(accepted_grid_points=len(accepted),singleton_grid_unresolved=singleton,
                            step_profile_bound_limited=row['design']=='A_step_step',
                            bracketing_grid_min=max(0.,min(accepted)-.025) if accepted else None,
                            bracketing_grid_max=min(.8,max(accepted)+.025) if accepted else None))
    summary['profile_resolution_and_boundary_flags']=flagged
    for row in summary['frozen_profiles']:
        relevant=[r for r in flagged if r['design']==row['design'] and r['sigma']==row['sigma']]
        row.update(singleton_grid_sets=sum(r['singleton_grid_unresolved'] for r in relevant),
                   width_interpretation='grid-unresolved, not zero uncertainty' if any(r['singleton_grid_unresolved'] for r in relevant) else
                   'upper endpoint imposed by fitted r bound; not d identification' if row['design']=='A_step_step' else 'grid extent only; threshold crossings not refined')
    summary['executed_source_snapshot']='inference_source_executed.py'
    summary['executed_source_sha256']=hashlib.sha256((HERE/'inference_source_executed.py').read_bytes()).hexdigest()
    summary['postprocessing_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    write_csv(HERE/'paired_sensitivity_summary.csv',paired)
    write_csv(HERE/'candidate_failures.csv',failure_rows)
    write_csv(HERE/'profile_resolution_flags.csv',flagged)
    write_csv(HERE/'parameter_recovery_summary.csv',summary['parameter_recovery'])
    write_csv(HERE/'frozen_profile_summary.csv',summary['frozen_profiles'])
    write_csv(HERE/'post_hoc_M1_paired_prediction.csv',m1_paired)
    (HERE/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(7.3,6.7),layout='constrained')
    confusion_source=[]
    for row, sigma in enumerate((.01,.03)):
        for col,design in enumerate(DESIGNS):
            ax=axes[row,col]
            matrix=np.array([[sum(r['generator']==g and r['candidate']==m for r in selected if r['sigma']==sigma and r['design']==design) for m in MODELS] for g in MODELS])
            ax.imshow(matrix,cmap='Blues',vmin=0,vmax=12)
            for i,g in enumerate(MODELS):
                for j,m in enumerate(MODELS):
                    ax.text(j,i,str(matrix[i,j]),ha='center',va='center',color='white' if matrix[i,j]>6 else 'black',fontsize=11)
                    confusion_source.append(dict(design=design,sigma=sigma,generator=g,selected=m,count=int(matrix[i,j])))
            ax.set(xticks=range(4),xticklabels=MODELS,yticks=range(4),yticklabels=MODELS,
                   xlabel='Training-selected family',ylabel='Generating family',title=f"{'A: step + step' if col==0 else 'B: step + combination'}; noise SD {sigma}")
    fig.suptitle('Family selection at equal measurement budgets',fontsize=12)
    finish(fig,'figure1_model_discrimination')
    write_csv(HERE/'figures/figure1_source.csv',confusion_source)
    fig,axes=plt.subplots(2,2,figsize=(9,6.6),layout='constrained')
    prediction_source=[]
    colors=('#0072B2','#D55E00')
    for row,sigma in enumerate((.01,.03)):
        for col,schedule in enumerate(('altered_pulse','full_withdrawal')):
            ax=axes[row,col]
            for gi,g in enumerate(MODELS):
                keys=sorted(set((r['case'],r['replicate']) for r in selected if r['sigma']==sigma and r['generator']==g))
                for ki,(case,rep) in enumerate(keys):
                    values=[]
                    jitter=(ki-(len(keys)-1)/2)*.008
                    for di,d in enumerate(DESIGNS):
                        r=next(v for v in selected if v['generator']==g and v['case']==case and v['replicate']==rep and v['sigma']==sigma and v['design']==d)
                        v=r[schedule+'_clean_rmse'];values.append(v)
                        prediction_source.append(dict(generator=g,case=case,sigma=sigma,replicate=rep,design=d,schedule=schedule,clean_rmse=v,selected=r['candidate']))
                    xx=np.array([gi-.17,gi+.17])+jitter
                    ax.plot(xx,values,color='#aaaaaa',alpha=.5,lw=.5,zorder=1)
                    ax.scatter(xx,values,c=colors,s=12,zorder=2)
                for di,d in enumerate(DESIGNS):
                    vv=[r[schedule+'_clean_rmse'] for r in selected if r['generator']==g and r['sigma']==sigma and r['design']==d]
                    ax.scatter([gi+(-.17 if di==0 else .17)],[np.median(vv)],marker='_',c=colors[di],s=170,linewidth=2.5,zorder=3)
            ax.axhline(sigma,color='#777777',ls=':',lw=1,label='Observation noise SD')
            ax.set(yscale='log',xticks=range(4),xticklabels=MODELS,ylabel='Clean held-out RMSE',xlabel='Generating family',title=f"{'Altered pulse' if col==0 else 'Full withdrawal'}; noise SD {sigma}")
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color=colors[0],label='A: step + step'),Line2D([0],[0],marker='o',ls='',color=colors[1],label='B: step + combination')],loc='outside upper center',ncol=2,frameon=False)
    finish(fig,'figure2_heldout_prediction')
    write_csv(HERE/'figures/figure2_source.csv',prediction_source)
    fig,axes=plt.subplots(1,3,figsize=(9,3.3),layout='constrained')
    profile_source=[]
    for i,d in enumerate(('A_step_step','B_step_combination','C_step_direct_d')):
        values=[v for v in profgrid if v['design']==d and float(v['sigma'])==.03 and int(v['replicate'])==0]
        x=[float(v['d']) for v in values];y=[max(0.,float(v['delta_chisq'])) for v in values]
        axes[i].plot(x,y,'o-',ms=3,color=('#777777','#D55E00','#009E73')[i])
        axes[i].axhline(3.841,color='black',ls='--',lw=.8)
        axes[i].axvline(.15,color='#777777',ls=':',lw=.8)
        axes[i].set(yscale='symlog',ylim=(-.1,1e5),xlim=(0,.8),xlabel='Fixed structural loss d',ylabel='Profile change in standardized SSE',title=('A: step + step','B: step + combination','C: step + direct d')[i])
        if i==0:axes[i].annotate('Upper limit from r bound',xy=(.8,y[-1]),xytext=(.21,200),arrowprops=dict(arrowstyle='->',lw=.8),fontsize=8)
        if i==2:axes[i].text(.02,8,'Narrow set unresolved\non 0.025 grid',fontsize=8)
        profile_source+=values
    fig.suptitle('Frozen M4; known a and b; noise SD 0.03, replicate 0',fontsize=11)
    finish(fig,'figure3_frozen_profiles')
    write_csv(HERE/'figures/figure3_source.csv',profile_source)
    print(json.dumps({'primary':summary['primary_dataset_count'],'candidate_failures':len(failure_rows),'true_family_failures':summary['true_family_failed_fit_count'],'posthoc':paired,'figures':3}))


if __name__=='__main__':main()
