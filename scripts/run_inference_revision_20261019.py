"""Prospective synthetic inference: noisy fits and whole-schedule holdouts.

No empirical calibration. CLI jobs are bounded and checkpoint candidate results.
Importing this module performs no scientific runs or file writes.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from scipy.integrate import odeint, ODEintWarning
from scipy.optimize import least_squares, minimize_scalar
from recoverability.model import Parameters, Schedule, solve

MODELS = ('M0', 'M1', 'M2', 'M3')
NAMES = {'M0': ('a', 'r'), 'M1': ('a', 'b', 'r'), 'M2': ('a', 'b', 'r'), 'M3': ('a', 'r', 'tau')}
RANGES = {'a': (.005, .2), 'b': (.001, 8.), 'r': (.2, 2.5), 'tau': (.2, 100.)}
TRAIN_T = np.unique(np.r_[np.arange(2., 121., 2.), .5, 1., 30.5, 31., 45., 45.5])
HOLD_T = np.unique(np.r_[np.arange(2., 121., 2.), .5, 1., 18., 27., 30.5, 31., 45., 45.5])
SCHEDULES = {
    'step': (30., 45., 0., None),
    'combination': (30., 45., .4, None),
    'altered_pulse': (18., 27., .25, None),
    'full_withdrawal': (30., 45., .4, 45.),
}
SIGMAS = (.01, .03)
DEFAULT_OUTPUT = ROOT / 'revision_20261019' / 'analysis'


def json_write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)


def cases():
    rows = np.random.default_rng(20261005).uniform([.03, 2., .8], [.09, 6., 1.3], (256, 3))[:2]
    return [dict(case='nominal', a=.06, b=4., r=1.05, tau=10.),
            *[dict(case=f'design_{i}', a=float(v[0]), b=float(v[1]), r=float(v[2]), tau=(1., 50.)[i]) for i, v in enumerate(rows)]]


def bounds(model):
    return np.log(np.array([RANGES[n] for n in NAMES[model]]).T)


def make_schedule(p, name):
    onset, end, amplitude, withdrawal = SCHEDULES[name]
    return Schedule(p['a'] + .08, p['a'], p['r'], onset=onset, pulse_end=end,
                    withdrawal=withdrawal, repair_amplitude=amplitude)


def core_reference(model, theta, schedule, times, frozen_d=None):
    names = NAMES[model]
    p = dict(zip(names, np.exp(theta)))
    par = Parameters(p['a'], p.get('b', 0.), p['r'], tau=p.get('tau', 10.))
    grid = np.r_[0., times]
    out = solve('M4' if frozen_d is not None else model, par, make_schedule(p, schedule),
                x0=.6, d0=0. if frozen_d is None else frozen_d, t_eval=grid)
    return 1. - out.x[1:]


def forward(model, theta, schedule, times, tolerance=1e-9, frozen_d=None):
    """y=1-x and dy/dlog(theta); all schedule boundaries are exact.

    LSODA independently duplicates the core equations. M0 has an analytic
    piecewise solution. Frozen M4 calls M1 with scaled repair and sensitivities.
    M3 uses its exact exponential input within each segment.
    """
    p = dict(zip(NAMES[model], np.exp(np.asarray(theta))))
    a, r, b, tau = p['a'], p['r'], p.get('b', 0.), p.get('tau', 10.)
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) == 0 or np.any(times <= 0) or np.any(np.diff(times) <= 0):
        raise ValueError('times must be strictly increasing and positive; x(0) is fixed')
    onset, end, amplitude, withdrawal = SCHEDULES[schedule]
    points = sorted(set([0., *[v for v in (onset, end, withdrawal) if v is not None and 0 < v < times[-1]], float(times[-1])]))
    k = len(theta)
    current = np.r_[.6, np.zeros(k)]
    out = np.empty((len(times), 1 + k))
    slow, slow_tau = a + .08, 0.
    scale = 1. if frozen_d is None else 1. - frozen_d
    for left, right in zip(points[:-1], points[1:]):
        midpoint = (left + right) / 2.
        active = midpoint >= onset and (withdrawal is None or midpoint < withdrawal)
        drive = a if active else a + .08
        repair = (r + (amplitude if active and midpoint < end else 0.)) * scale
        indices = np.flatnonzero((times > left) & (times <= right))
        grid = np.unique(np.r_[left, times[indices], right])
        if model == 'M0' and frozen_d is None:
            dt = grid - left
            rate = drive + repair
            exp = np.exp(-rate*dt)
            eq = drive / rate
            drate = np.array([a, r])
            deq = np.array([a*repair, -drive*r]) / rate**2
            states = np.empty((len(grid), 3))
            states[:, 0] = eq + (current[0]-eq)*exp
            states[:, 1:] = deq + (current[1:]-deq)*exp[:, None] - (current[0]-eq)*exp[:, None]*dt[:, None]*drate
        else:
            def derivative(z, tt):
                x = z[0]
                if model == 'M3':
                    dt = tt - left
                    e = math.exp(-dt/tau)
                    s = drive + (slow-drive)*e
                    st = slow_tau*e + (slow-drive)*e*dt/tau
                    f = s*(1.-x)-repair*x
                    fx = -s-repair
                    forcing = (a*(1.-x), -r*x, st*(1.-x))
                elif model == 'M2':
                    f = drive*(1.-x)+b*x*(1.-x)-repair*x
                    fx = -drive+b*(1.-2.*x)-repair
                    forcing = (a*(1.-x), b*x*(1.-x), -r*scale*x)
                else:
                    f = drive*(1.-x)+b*x*x*(1.-x)-repair*x
                    fx = -drive+b*(2.*x-3.*x*x)-repair
                    forcing = (a*(1.-x), b*x*x*(1.-x), -r*scale*x)
                return [f, *[fx*z[j+1]+forcing[j] for j in range(k)]]
            with warnings.catch_warnings():
                warnings.simplefilter('error', ODEintWarning)
                states, info = odeint(derivative, current, grid, rtol=tolerance, atol=tolerance*.01,
                                     mxstep=5000, full_output=True)
            if info['message'] != 'Integration successful.':
                raise ArithmeticError(info['message'])
        if not np.isfinite(states).all() or np.min(states[:, 0]) < -1e-8 or np.max(states[:, 0]) > 1.+1e-8:
            raise ArithmeticError('nonfinite/unphysical latent state; no clipping')
        if len(indices):
            out[indices] = states[np.searchsorted(grid, times[indices])]
        current = states[-1].copy()
        if model == 'M3':
            duration = right-left
            e = math.exp(-duration/tau)
            slow_tau = slow_tau*e+(slow-drive)*e*duration/tau
            slow = drive+(slow-drive)*e
    return 1.-out[:, 0], -out[:, 1:]


def fit_candidate(model, training, sigma):
    low, high = bounds(model)
    starts, best = [], None
    for start_index, fraction in enumerate((.2, .5, .8)):
        t0 = time.monotonic()
        cache = {'theta': None}
        def evaluate(theta):
            if cache['theta'] is None or not np.array_equal(theta, cache['theta']):
                packets = [forward(model, theta, schedule, TRAIN_T) for schedule, _ in training]
                cache.update(theta=theta.copy(), residual=np.concatenate([(out[0]-obs)/sigma for out, (_, obs) in zip(packets, training)]),
                             jac=np.vstack([out[1]/sigma for out in packets]))
            return cache
        try:
            fit = least_squares(lambda q: evaluate(q)['residual'], low+fraction*(high-low),
                                jac=lambda q: evaluate(q)['jac'], bounds=(low, high), max_nfev=150,
                                ftol=1e-8, xtol=1e-8, gtol=1e-6)
            rss = float(np.dot(fit.fun, fit.fun))
            entry = dict(start=start_index, start_fraction=fraction, success=bool(fit.success), status=int(fit.status),
                         message=str(fit.message), nfev=int(fit.nfev), optimality=float(fit.optimality),
                         standardized_sse=rss, theta_log=fit.x.tolist(), seconds=time.monotonic()-t0)
            if fit.success and (best is None or rss < best['standardized_sse']):
                best = entry.copy()
        except Exception as exc:
            entry = dict(start=start_index, start_fraction=fraction, success=False, message=f'{type(exc).__name__}: {exc}', seconds=time.monotonic()-t0)
        starts.append(entry)
    result = dict(candidate=model, success=best is not None, starts=starts)
    if best is not None:
        theta = np.array(best['theta_log'])
        result.update(best)
        result.update(parameters=dict(zip(NAMES[model], np.exp(theta).tolist())),
                      bic=best['standardized_sse']+len(theta)*math.log(len(TRAIN_T)*len(training)),
                      boundary_parameters=[name for name, q, lo, hi in zip(NAMES[model], theta, low, high) if min(q-lo, hi-q) <= .01*(hi-lo)],
                      start_loss_spread=max(s['standardized_sse'] for s in starts if s['success'])-min(s['standardized_sse'] for s in starts if s['success']))
    return result


def noise(seed, model_index, case_index, sigma_index, replicate, trace, size):
    rng = np.random.default_rng(np.random.SeedSequence([seed, model_index, case_index, sigma_index, replicate, trace]))
    return rng.normal(0., SIGMAS[sigma_index], size)


def provenance(mode):
    return dict(evidence_class='synthetic_simulation', empirical_calibration=False, mode=mode,
                started_utc=datetime.now(timezone.utc).isoformat(),
                source_hashes={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), ROOT/'protocol/inference_revision_20261019.md', ROOT/'src/recoverability/model.py']},
                versions={n: importlib.metadata.version(n) for n in ('numpy', 'scipy')}, python=platform.python_version(),
                parameters=cases(), training_times=TRAIN_T.tolist(), heldout_times=HOLD_T.tolist(), schedules=SCHEDULES,
                x0_known=.6, optimization_starts=[.2, .5, .8], known_noise_sigmas=list(SIGMAS),
                n_per_training_design=len(TRAIN_T)*2, replicates=4, status='running')


def validate(output):
    rows = []
    for model in MODELS:
        low, high = bounds(model)
        settings = [np.log([c[n] for n in NAMES[model]]) for c in cases()]
        settings += [low+.1*(high-low), low+.9*(high-low)]
        for setting_index, theta in enumerate(settings):
            for schedule in SCHEDULES:
                y, jac = forward(model, theta, schedule, HOLD_T)
                ref = core_reference(model, theta, schedule, HOLD_T)
                fine, _ = forward(model, theta, schedule, HOLD_T, tolerance=1e-11)
                fd = np.column_stack([(forward(model, theta+np.eye(len(theta))[j]*1e-5, schedule, HOLD_T)[0]-
                                       forward(model, theta-np.eye(len(theta))[j]*1e-5, schedule, HOLD_T)[0])/2e-5 for j in range(len(theta))])
                err = float(np.max(np.abs(y-ref)))
                refine = float(np.max(np.abs(y-fine)))
                je = float(np.max(np.abs(jac-fd)/(1.+np.abs(fd))))
                rows.append(dict(model=model, setting=setting_index, schedule=schedule, max_reference_error=err, max_refinement_error=refine, max_jacobian_scaled_error=je,
                                 passed=err<=1e-6 and refine<=1e-6 and je<=1e-4))
    frozen_theta = np.log([.06, 4., 1.05])
    for d in (0., .15, .8):
        for schedule in SCHEDULES:
            y, _ = forward('M1', frozen_theta, schedule, HOLD_T, frozen_d=d)
            ref = core_reference('M1', frozen_theta, schedule, HOLD_T, frozen_d=d)
            err = float(np.max(np.abs(y-ref)))
            rows.append(dict(model='frozen_M4', setting=d, schedule=schedule, max_reference_error=err, passed=err<=1e-6))
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output/'numerical_validation.csv', rows)
    result = dict(passed=all(r['passed'] for r in rows), checks=len(rows), max_reference_error=max(r['max_reference_error'] for r in rows),
                  max_jacobian_scaled_error=max(r.get('max_jacobian_scaled_error',0.) for r in rows), provenance=provenance('validation'))
    json_write(output/'numerical_validation.json', result)
    return result


def run_generator(model, output, timeout=600, smoke=False):
    validation = json.loads((output/'numerical_validation.json').read_text('utf-8'))
    if not validation['passed']:
        raise ArithmeticError('Numerical validation must pass before fits')
    directory = output/('smoke' if smoke else model)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory/'manifest.json').exists() and not smoke:
        raise FileExistsError('Refusing to replace an existing job; preserve complete/failed scientific records')
    manifest = provenance('smoke' if smoke else model)
    json_write(directory/'manifest.json', manifest)
    rows, observation_rows, predictions, failures = [], [], [], []
    model_index = MODELS.index(model)
    started = time.monotonic()
    interrupted = False
    for ci, case in enumerate(cases()[:1] if smoke else cases()):
        true_theta = np.log([case[n] for n in NAMES[model]])
        clean_training = {s: forward(model,true_theta,s,TRAIN_T)[0] for s in ('step','combination')}
        clean_hold = {s: forward(model,true_theta,s,HOLD_T)[0] for s in ('altered_pulse','full_withdrawal')}
        for si, sigma in enumerate(SIGMAS[:1] if smoke else SIGMAS):
            for rep in range(1 if smoke else 4):
                first_error = noise(20261019,model_index,ci,si,rep,0,len(TRAIN_T))
                second_error = noise(20261019,model_index,ci,si,rep,1,len(TRAIN_T))
                held_noisy = {s: value+noise(20261020,model_index,ci,si,rep,j,len(HOLD_T)) for j,(s,value) in enumerate(clean_hold.items())}
                for design, second_schedule in (('A_step_step','step'),('B_step_combination','combination')):
                    if time.monotonic()-started >= timeout-15:
                        failures.append(dict(reason='bounded job runtime reached',case=case['case'],sigma=sigma,replicate=rep,design=design)); interrupted=True; break
                    key = dict(generator=model,case=case['case'],sigma=sigma,replicate=rep,design=design)
                    training = [('step',clean_training['step']+first_error),(second_schedule,clean_training[second_schedule]+second_error)]
                    for trace,(schedule,observed) in enumerate(training):
                        for t, clean, noisy in zip(TRAIN_T,clean_training[schedule],observed):
                            observation_rows.append(dict(**key,role='training',trace=trace,schedule=schedule,time=float(t),clean=float(clean),observed=float(noisy)))
                    for h,schedule in enumerate(clean_hold):
                        for t,clean,noisy in zip(HOLD_T,clean_hold[schedule],held_noisy[schedule]):
                            observation_rows.append(dict(**key,role='heldout',trace=h,schedule=schedule,time=float(t),clean=float(clean),observed=float(noisy)))
                    fitted = [fit_candidate(candidate,training,sigma) for candidate in MODELS]
                    successful = sorted([f for f in fitted if f['success']], key=lambda f:f['bic'])
                    winner = successful[0]['candidate'] if successful else None
                    gap = successful[1]['bic']-successful[0]['bic'] if len(successful)>1 else None
                    for result in fitted:
                        record = {**key,**result,'selected':result['candidate']==winner,'runner_up_bic_gap':gap,'ambiguous_selection':gap is not None and gap<2.}
                        if result['success']:
                            theta = np.array(result['theta_log'])
                            numerical_pass = True
                            for schedule in clean_hold:
                                predicted,_ = forward(result['candidate'],theta,schedule,HOLD_T)
                                reference = core_reference(result['candidate'],theta,schedule,HOLD_T)
                                err = float(np.max(np.abs(predicted-reference)))
                                numerical_pass &= err<=1e-6
                                record[schedule+'_reference_error'] = err
                                record[schedule+'_clean_rmse'] = float(np.sqrt(np.mean((predicted-clean_hold[schedule])**2)))
                                record[schedule+'_noisy_rmse'] = float(np.sqrt(np.mean((predicted-held_noisy[schedule])**2)))
                                for t,clean,obs,pred in zip(HOLD_T,clean_hold[schedule],held_noisy[schedule],predicted):
                                    predictions.append(dict(**key,candidate=result['candidate'],selected=record['selected'],schedule=schedule,time=float(t),clean=float(clean),observed=float(obs),prediction=float(pred)))
                            record['numerical_pass'] = bool(numerical_pass)
                            if result['candidate']==model:
                                record['parameter_relative_errors'] = {n:abs(v-case[n])/case[n] for n,v in result['parameters'].items()}
                                # Diagnostic only, calculated after all candidate fits and selection.
                                # The feasible truth is never an initialization or replacement fit.
                                truth_sse = float(sum(np.sum(((clean_training[s]-obs)/sigma)**2) for s,obs in training))
                                record['diagnostic_truth_standardized_sse'] = truth_sse
                                record['diagnostic_fit_minus_truth_sse'] = result['standardized_sse']-truth_sse
                                record['diagnostic_optimizer_miss'] = bool(result['standardized_sse']>truth_sse+1e-3)
                        rows.append(record)
                        with (directory/'fit_checkpoints.jsonl').open('a',encoding='utf-8') as f:
                            f.write(json.dumps(record,allow_nan=False)+'\n')
                    print(json.dumps(dict(progress=key,winner=winner,gap=gap,elapsed_seconds=round(time.monotonic()-started,1))),flush=True)
                if interrupted: break
            if interrupted: break
        if interrupted: break
    json_write(directory/'fits.json', rows)
    flat = []
    for r in rows:
        q = {k:v for k,v in r.items() if not isinstance(v,(dict,list))}
        q.update({'estimate_'+k:v for k,v in r.get('parameters',{}).items()})
        q.update({'relative_error_'+k:v for k,v in r.get('parameter_relative_errors',{}).items()})
        q['boundary_parameters'] = ';'.join(r.get('boundary_parameters',[]))
        flat.append(q)
    write_csv(directory/'fits.csv',flat)
    write_csv(directory/'observations.csv',observation_rows)
    write_csv(directory/'heldout_predictions.csv',predictions)
    expected = 8 if smoke else 192
    manifest.update(status='complete_unreviewed' if len(rows)==expected and not failures else 'incomplete_unreviewed',
                    finished_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started,
                    expected_candidate_fits=expected,actual_candidate_fits=len(rows),failures=failures,
                    file_hashes={p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!='manifest.json'})
    json_write(directory/'manifest.json',manifest)
    return manifest


def frozen_prediction(r,d,schedule):
    return forward('M1',np.log([.06,4.,r]),schedule,TRAIN_T,frozen_d=d)[0]


def run_profiles(output,timeout=600):
    directory=output/'frozen_M4'
    directory.mkdir(parents=True,exist_ok=True)
    if (directory/'manifest.json').exists():
        raise FileExistsError('Refusing to overwrite frozen profile job')
    manifest=provenance('frozen_M4_profiles')
    manifest.update(known_parameters=dict(a=.06,b=4.,epsilon=0.,eta=0.,x0=.6),truth=dict(r=1.05,d=.15),profile_grid=np.linspace(0.,.8,33).tolist())
    json_write(directory/'manifest.json',manifest)
    started=time.monotonic(); rows=[]; fitted=[]; observations=[]
    truth_step=frozen_prediction(1.05,.15,'step'); truth_combo=frozen_prediction(1.05,.15,'combination')
    for si,sigma in enumerate(SIGMAS):
        for rep in range(4):
            rng=np.random.default_rng(np.random.SeedSequence([20261021,si,rep]))
            errors=rng.normal(0.,sigma,(2,len(TRAIN_T)))
            for design in ('A_step_step','B_step_combination','C_step_direct_d'):
                if time.monotonic()-started>timeout-15:
                    raise TimeoutError('bounded profile runtime reached; prior profile checkpoints retained')
                y1=truth_step+errors[0]
                y2=(truth_step if design=='A_step_step' else truth_combo if design=='B_step_combination' else np.full(len(TRAIN_T),.15))+errors[1]
                def residual(r,d):
                    p1=frozen_prediction(r,d,'step')
                    p2=p1 if design=='A_step_step' else frozen_prediction(r,d,'combination') if design=='B_step_combination' else np.full(len(TRAIN_T),d)
                    return np.r_[p1-y1,p2-y2]/sigma
                start_results=[]
                for rr,dd in ((.4,.1),(1.5,.4),(3.,.7)):
                    fit=least_squares(lambda q:residual(math.exp(q[0]),q[1]),[math.log(rr),dd],bounds=([math.log(.2),0.],[math.log(4.),.8]),max_nfev=150,ftol=1e-8,xtol=1e-8,gtol=1e-6)
                    start_results.append(dict(success=bool(fit.success),r=math.exp(fit.x[0]),d=float(fit.x[1]),sse=float(fit.fun@fit.fun),nfev=int(fit.nfev),optimality=float(fit.optimality)))
                valid=[v for v in start_results if v['success']]
                best=min(valid,key=lambda v:v['sse']) if valid else None
                key=dict(design=design,sigma=sigma,replicate=rep)
                local=[]
                for d in np.linspace(0.,.8,33):
                    opt=minimize_scalar(lambda r:float(residual(r,d)@residual(r,d)),bounds=(.2,4.),method='bounded',options=dict(xatol=1e-9,maxiter=150))
                    local.append(dict(**key,d=float(d),fitted_r=float(opt.x),sse=float(opt.fun),success=bool(opt.success),nfev=int(opt.nfev)))
                minimum=min([v['sse'] for v in local if v['success']]+([best['sse']] if best else []))
                for v in local: v.update(delta_chisq=v['sse']-minimum,in_profile_set=v['sse']-minimum<=3.841)
                interval=[v['d'] for v in local if v['in_profile_set']]
                row=dict(**key,starts=start_results,success=best is not None,profile_set_min=min(interval) if interval else None,profile_set_max=max(interval) if interval else None,
                         profile_set_width=max(interval)-min(interval) if interval else None,profile_grid_size=len(local),profile_sse_range=max(v['sse'] for v in local)-min(v['sse'] for v in local))
                if best: row.update(best,r_relative_error=abs(best['r']-1.05)/1.05,d_absolute_error=abs(best['d']-.15))
                rows+=local;fitted.append(row)
                for j,(clean,obs,channel) in enumerate(((truth_step,y1,'anchored_1_minus_x'),(truth_step if design=='A_step_step' else truth_combo if design=='B_step_combination' else np.full(len(TRAIN_T),.15),y2,'direct_d' if design=='C_step_direct_d' else 'anchored_1_minus_x'))):
                    observations.extend(dict(**key,trace=j,channel=channel,time=float(t),clean=float(c),observed=float(o)) for t,c,o in zip(TRAIN_T,clean,obs))
                with (directory/'checkpoints.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(dict(fit=row,profile=local))+'\n')
                print(json.dumps(dict(progress=key,profile_width=row['profile_set_width'],elapsed_seconds=round(time.monotonic()-started,1))),flush=True)
    write_csv(directory/'profiles.csv',rows);write_csv(directory/'observations.csv',observations);json_write(directory/'fits.json',fitted)
    write_csv(directory/'fits.csv',[{k:v for k,v in r.items() if k!='starts'} for r in fitted])
    manifest.update(status='complete_unreviewed',elapsed_seconds=time.monotonic()-started,actual_fit_datasets=len(fitted),actual_profile_points=len(rows),finished_utc=datetime.now(timezone.utc).isoformat(),
                    file_hashes={p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!='manifest.json'})
    json_write(directory/'manifest.json',manifest)
    return manifest


def summarize(output):
    fits=[]; manifests=[]
    for model in MODELS:
        directory=output/model
        fits+=json.loads((directory/'fits.json').read_text('utf-8'))
        manifests.append(json.loads((directory/'manifest.json').read_text('utf-8')))
    summaries=[]
    for design in ('A_step_step','B_step_combination'):
        for sigma in SIGMAS:
            for generator in (*MODELS,'all'):
                selected=[r for r in fits if r['design']==design and r['sigma']==sigma and r['selected'] and (generator=='all' or r['generator']==generator)]
                correct=[r for r in selected if r['candidate']==r['generator']]
                row=dict(design=design,sigma=sigma,generator=generator,n=len(selected),correct=len(correct),ambiguous=sum(r['ambiguous_selection'] for r in selected),
                         selected_numerical_failures=sum(not r['numerical_pass'] for r in selected),confusion={m:sum(r['candidate']==m for r in selected) for m in MODELS})
                for schedule in ('altered_pulse','full_withdrawal'):
                    for kind in ('clean','noisy'):
                        vals=[r[f'{schedule}_{kind}_rmse'] for r in selected]
                        row[f'{schedule}_{kind}_rmse_median']=float(np.median(vals))
                        row[f'{schedule}_{kind}_rmse_max']=float(max(vals))
                summaries.append(row)
    param=[]
    for model in MODELS:
        for design in ('A_step_step','B_step_combination'):
            for sigma in SIGMAS:
                true=[r for r in fits if r['generator']==model and r['candidate']==model and r['design']==design and r['sigma']==sigma and r['success']]
                for name in NAMES[model]:
                    values=[r['parameter_relative_errors'][name] for r in true]
                    param.append(dict(generator=model,design=design,sigma=sigma,parameter=name,n=len(true),median_relative_error=float(np.median(values)),max_relative_error=float(max(values)),boundary_fits=sum(name in r['boundary_parameters'] for r in true)))
    profilefits=json.loads((output/'frozen_M4/fits.json').read_text('utf-8'))
    profilesummary=[]
    for design in ('A_step_step','B_step_combination','C_step_direct_d'):
        for sigma in SIGMAS:
            rows=[r for r in profilefits if r['design']==design and r['sigma']==sigma]
            profilesummary.append(dict(design=design,sigma=sigma,n=len(rows),median_profile_width=float(np.median([r['profile_set_width'] for r in rows])),
                                       median_r_relative_error=float(np.median([r['r_relative_error'] for r in rows])),median_d_absolute_error=float(np.median([r['d_absolute_error'] for r in rows])),
                                       median_profile_sse_range=float(np.median([r['profile_sse_range'] for r in rows]))))
    starts=[s for r in fits for s in r['starts']]
    result=dict(evidence_class='synthetic_simulation',empirical_calibration=False,status='generated_pending_independent_review',
                primary_candidate_fits=len(fits),primary_dataset_count=len(fits)//4,optimizer_starts=len(starts),failed_starts=sum(not s['success'] for s in starts),
                failed_candidate_fits=sum(not r['success'] for r in fits),candidate_numerical_failures=sum(not r.get('numerical_pass',False) for r in fits if r['success']),
                maximum_fitted_holdout_reference_error=max(r.get(s+'_reference_error',0.) for r in fits for s in ('altered_pulse','full_withdrawal')),
                any_boundary_candidate_fits=sum(bool(r.get('boundary_parameters')) for r in fits),selections=summaries,parameter_recovery=param,frozen_profiles=profilesummary,
                diagnostic_true_family_optimizer_misses=sum(r.get('diagnostic_optimizer_miss',False) for r in fits),
                primary_elapsed_seconds=sum(m['elapsed_seconds'] for m in manifests))
    json_write(output/'summary.json',result)
    write_csv(output/'selection_summary.csv',[{k:(json.dumps(v) if isinstance(v,dict) else v) for k,v in r.items()} for r in summaries])
    write_csv(output/'parameter_recovery_summary.csv',param);write_csv(output/'frozen_profile_summary.csv',profilesummary)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('validate','smoke','run','profiles','summarize'),required=True)
    parser.add_argument('--generator',choices=MODELS,default='M1')
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--timeout',type=int,default=600)
    args=parser.parse_args()
    if not 1<=args.timeout<=600:parser.error('timeout must be1..600seconds')
    if args.mode=='validate':result=validate(args.output)
    elif args.mode=='smoke':result=run_generator(args.generator,args.output,args.timeout,True)
    elif args.mode=='run':result=run_generator(args.generator,args.output,args.timeout)
    elif args.mode=='profiles':result=run_profiles(args.output,args.timeout)
    else:result=summarize(args.output)
    print(json.dumps(result,allow_nan=False))
    return 0 if result.get('passed',True) and result.get('status')!='incomplete_unreviewed' else 1


if __name__=='__main__':
    raise SystemExit(main())
