"""Bounded synthetic study driver. --help/import have no scientific side effects."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time
import traceback

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261005
THRESHOLDS = (.10, .15, .20)
SAMPLE_STEP = .25
MAX_STEP = .5
ARMS = (
    ("none", "identical_policies", False, 0., None),
    ("repair", "identical_policies", False, .4, None),
    ("drive", "maintenance", True, 0., None),
    ("drive", "full_withdrawal", True, 0., 45.),
    ("both", "maintenance", True, .4, None),
    ("both", "full_withdrawal", True, .4, 45.),
)
MODELS = (("M0", 10.), ("M1", 10.), ("M2", 10.),
          ("M3", 1.), ("M3", 10.), ("M3", 50.))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csv_write(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


class Budget:
    def __init__(self, seconds):
        self.start = time.monotonic()
        self.seconds = seconds

    def check(self):
        if time.monotonic() - self.start > self.seconds:
            raise TimeoutError(f"Job budget exceeded ({self.seconds} seconds)")


def roots(model, a, b, r, d=0.):
    """Independent polynomial reference, explicitly retaining degeneracies."""
    re = r * (1-d) if model == "M4" else r
    coefficients = ([-b, b, -a-re, a] if model in ("M1", "M4") else
                    [-b, b-a-re, a] if model == "M2" else [-a-re, a])
    coefficients = np.trim_zeros(np.asarray(coefficients, float), "f")
    if not len(coefficients):
        return [], True
    found = []
    for z in np.roots(coefficients):
        if abs(z.imag) > 1e-7 or not -1e-9 <= z.real <= 1+1e-9:
            continue
        x = float(z.real)
        residual = abs(float(np.polyval(coefficients, x)))
        slope = float(np.polyval(np.polyder(coefficients), x))
        kind = "nonhyperbolic" if abs(slope) <= 1e-7 else "stable" if slope < 0 else "unstable"
        if residual >= 1e-10:
            kind = "uncertified"
        if found and any(abs(x-p["x"]) < 1e-6 for p in found):
            for p in found:
                if abs(x-p["x"]) < 1e-6:
                    p["stability"] = "nonhyperbolic"
            continue
        found.append(dict(x=x, derivative=slope, stability=kind, residual=residual))
    return sorted(found, key=lambda p: p["x"]), False


def setting(a=.06, b=4., r=1.05, **extra):
    return dict(a=a, b=b, r=r, d0=0., epsilon=0., eta=0., xc=.25, **extra)


def initial_x(model, p):
    eq, continuum = roots(model, p["a"]+.08, p["b"], p["r"], p.get("d0", 0.))
    stable = [q["x"] for q in eq if q["stability"] == "stable"]
    if continuum or not stable:
        raise ValueError("No residual-certified attracting pretreatment root")
    return max(stable)


def design(quick=False):
    rng = np.random.default_rng(SEED)
    draws = rng.uniform([.03, 2., .8], [.09, 6., 1.3], (256, 3))
    rows = [setting(setting_id="nominal", design_kind="nominal")]
    rows.extend(setting(*map(float, values), setting_id=f"sample_{i:03}", design_kind="uniform_design")
                for i, values in enumerate(draws))
    return rows[:3] if quick else rows


def make_cases(model, tau, parameters, x_override=None):
    cases = []
    for p in parameters:
        for arm, policy, reduce, amplitude, withdrawal in ARMS:
            cases.append(dict(p, model=model, tau=tau, arm=arm, policy=policy,
                              reduce_drive=reduce, amplitude=amplitude, withdrawal=withdrawal,
                              pulse_end=45., x0=initial_x(model, p) if x_override is None else x_override,
                              initial_selection="highest_attracting_pretreatment_root" if x_override is None else "prespecified_x0",
                              b_used=model not in ("M0", "M3"), initial_state_reference="t0"))
    return cases


def vector_rhs(model, a, b, r, tau, epsilon, eta, xc, y):
    x, s, d = y[:, 0], y[:, 1], y[:, 2]
    dx = (s if model == "M3" else a) * (1-x) - r*x*(1-d if model == "M4" else 1)
    if model in ("M1", "M4"):
        dx = dx + b*x*x*(1-x)
    elif model == "M2":
        dx = dx + b*x*(1-x)
    ds = (a-s)/tau if model == "M3" else np.zeros_like(x)
    dd = epsilon*np.maximum(x-xc, 0)*(1-d)-eta*d if model == "M4" else np.zeros_like(x)
    return np.column_stack((dx, ds, dd))


def integrate(model, cases, horizon, budget, ramp=None):
    """Independent vectorized DOP853, same equations and exact event segments."""
    n = len(cases)
    array = lambda key: np.array([p[key] for p in cases], float)
    a0, b, r0, tau = array("a"), array("b"), array("r"), array("tau")
    eps, eta, xc = array("epsilon"), array("eta"), array("xc")
    pre = a0+.08
    repair, ends = array("amplitude"), array("pulse_end")
    reduce = np.array([p["reduce_drive"] for p in cases])
    withdrawal = np.array([np.inf if p["withdrawal"] is None else p["withdrawal"] for p in cases])
    start_time = 30. if all(p.get("initial_state_reference")=="treatment_onset" for p in cases) else 0.
    initial_s = np.array([p.get("s0",p.get("fixed_command",p["a"]+.08)) for p in cases],float)
    y = np.column_stack((array("x0"), initial_s, array("d0")))
    if ramp:
        y[:, 1] = .02
    times = np.unique(np.r_[np.arange(0, horizon+SAMPLE_STEP/2, SAMPLE_STEP), [30., 45., horizon], ends, withdrawal[np.isfinite(withdrawal)]])
    times = times[(times >= start_time) & (times <= horizon)]
    events = sorted(set([start_time, horizon] + ([] if ramp or start_time>=30 else [30.]) +
                        [float(v) for v in np.r_[ends, withdrawal] if start_time < v < horizon]))
    if ramp:
        events = [0., float(ramp), 2.*float(ramp)]
        times = np.unique(np.r_[np.linspace(0, horizon, min(20001, int(horizon/.25)+1)), events])
    result = np.full((len(times), n, 3), np.nan)
    result[0] = y
    for left, right in zip(events[:-1], events[1:]):
        budget.check()
        middle = (left+right)/2
        active_drive = reduce & (middle >= 30) & (middle < withdrawal)
        aa = np.array([p.get("fixed_command",v) for p,v in zip(cases,np.where(active_drive, a0, pre))],float)
        rr = r0 + np.where((middle >= 30) & (middle < ends) & (middle < withdrawal), repair, 0.)
        def fun(t, flat):
            budget.check()
            command = aa if ramp is None else np.full(n, .02+.12*(t/ramp if t <= ramp else 2-t/ramp))
            return vector_rhs(model, command, b, rr if ramp is None else r0, tau, eps, eta, xc, flat.reshape(n, 3)).ravel()
        sol = solve_ivp(fun, (left, right), y.ravel(), method="DOP853", rtol=1e-9, atol=1e-11, max_step=MAX_STEP, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        internal=sol.y.T.reshape(-1,n,3)
        if not np.isfinite(internal).all() or internal[:,:,1].min() < -1e-9:
            raise RuntimeError("Internal solver states nonfinite or input state negative")
        if internal[:,:,[0,2]].min() < -1e-9 or internal[:,:,[0,2]].max() > 1+1e-9:
            raise RuntimeError("Internal solver state bounds exceeded; no clipping permitted")
        mask = (times > left) & (times <= right)
        result[mask] = sol.sol(times[mask]).T.reshape(-1, n, 3)
        y = sol.y[:, -1].reshape(n, 3)
    if not np.isfinite(result).all():
        raise RuntimeError("Nonfinite/unfilled trajectory output")
    if result[:,:,1].min() < -1e-9:
        raise RuntimeError("Sampled input state negative")
    if result[:, :, [0, 2]].min() < -1e-9 or result[:, :, [0, 2]].max() > 1+1e-9:
        raise RuntimeError("State bounds exceeded; no clipping permitted")
    return times, result


def outcome_rows(model, cases, t, y, horizon):
    rows = []
    post = t >= 30
    persist = (t >= horizon-30) & (t <= horizon)
    last = t >= horizon-10
    for j, p in enumerate(cases):
        command = p.get("fixed_command",p["a"] if p["reduce_drive"] and p["withdrawal"] is None else p["a"]+.08)
        drift = vector_rhs(model, command, p["b"], p["r"], p["tau"], p["epsilon"], p["eta"], p["xc"], y[last, j])
        max_drift = float(np.max(np.abs(drift)))
        eq, continuum = roots(model, command, p["b"], p["r"], y[-1, j, 2])
        attracting = [q["x"] for q in eq if q["stability"] == "stable"]
        distance = min((abs(y[-1, j, 0]-v) for v in attracting), default=np.inf)
        if model == "M4" and (p["epsilon"] or p["eta"]):
            convergence = "full_drift_small_equilibrium_certification_deferred" if max_drift < 1e-8 else "unresolved_at_horizon"
        else:
            convergence = "converged_frozen_structure_manifold" if model == "M4" and max_drift < 1e-8 and distance < 1e-6 else "converged" if max_drift < 1e-8 and distance < 1e-6 else "unresolved_at_horizon"
        actual_cessation = min(p["pulse_end"],p["withdrawal"] if p["withdrawal"] is not None else np.inf)
        cessation = int(np.argmin(abs(t-actual_cessation)))
        ceq, _ = roots(model, command, p["b"], p["r"], y[cessation, j, 2])
        unstable = [q["x"] for q in ceq if q["stability"] == "unstable"]
        nearest_unstable = min(unstable, key=lambda v: abs(v-y[cessation, j, 0])) if unstable else np.nan
        xmax = float(y[persist, j, 0].max())
        dmax = float(y[persist, j, 2].max())
        prepoints,_=roots(model,p.get("fixed_command",p["a"]+.08),p["b"],p["r"],p["d0"])
        postpoints,_=roots(model,p.get("fixed_command",p["a"]),p["b"],p["r"],p["d0"])
        pre_stable=[q["x"] for q in prepoints if q["stability"]=="stable"]
        post_stable=[q["x"] for q in postpoints if q["stability"]=="stable"]
        row = dict(p, endpoint=horizon, x_endpoint=float(y[-1, j, 0]), s_endpoint=float(y[-1, j, 1]), d_endpoint=float(y[-1, j, 2]),
                   q_endpoint=float((1-y[-1, j, 0])*(1-y[-1, j, 2])), x_persistence_max=xmax, d_persistence_max=dmax,
                   impairment_auc=float(np.trapezoid(y[post, j, 0], t[post])), max_drift_last10=max_drift,
                   attracting_root_distance=float(distance), convergence=convergence,
                   x_at_cessation=float(y[cessation, j, 0]), unstable_root_at_cessation=nearest_unstable,
                   signed_distance_to_unstable=float(y[cessation, j, 0]-nearest_unstable),
                   persistence_evaluation="sampled_every_.25_plus_events", root_continuum=continuum,
                   actual_repair_cessation=actual_cessation,actual_repair_duration=max(0.,actual_cessation-30.),
                   x30=float(y[np.argmin(abs(t-30)),j,0]), pre_root_count=len(prepoints),post_root_count=len(postpoints),
                   pre_stabilities="|".join(q["stability"] for q in prepoints),post_stabilities="|".join(q["stability"] for q in postpoints),
                   pre_bistable=[q["stability"] for q in prepoints]==["stable","unstable","stable"],
                   post_bistable=[q["stability"] for q in postpoints]==["stable","unstable","stable"],
                   frozen_pre_target_root_exists=any(x<=.15 for x in pre_stable),
                   frozen_post_target_root_exists=any(x<=.15 for x in post_stable),
                   equilibrium_flag_scope="frozen_d_only" if model=="M4" else "settled_input" if model=="M3" else "scalar")
        for threshold in THRESHOLDS:
            tag = f"{threshold:.2f}"
            row[f"x_target_{tag}"] = xmax <= threshold
            row[f"d_target_{tag}"] = dmax <= threshold
            row[f"joint_target_{tag}"] = xmax <= threshold and (model != "M4" or dmax <= threshold)
            row[f"x_target_margin_{tag}"] = threshold-xmax
            row[f"d_target_margin_{tag}"] = threshold-dmax
            row[f"x_target_borderline_{tag}"] = abs(threshold-xmax)<=1e-6
            row[f"d_target_borderline_{tag}"] = model=="M4" and abs(threshold-dmax)<=1e-6
            row[f"frozen_pre_target_root_exists_{tag}"] = any(x<=threshold for x in pre_stable)
            row[f"frozen_post_target_root_exists_{tag}"] = any(x<=threshold for x in post_stable)
        row["primary_target"] = row["joint_target_0.15"]
        drive_end = p["withdrawal"] if p["withdrawal"] is not None else horizon
        row["cost_drive"] = max(0., min(horizon, drive_end)-30) if p["reduce_drive"] else 0.
        row["cost_repair"] = max(0., min(horizon, actual_cessation)-30)*(p["amplitude"]/.4)**2
        row["cost_total"] = row["cost_drive"]+row["cost_repair"]
        if model == "M4":
            row["baseline_d_meets_target_only"] = p["d0"] <= .15
            row["baseline_flag_not_reachability"] = True
            row["structural_repair_possible"] = p["eta"] > 0
        rows.append(row)
    return rows


def save_nominal(out, stem, cases, t, y):
    selected = [j for j,p in enumerate(cases) if p["setting_id"] == "nominal"]
    csv_write(out/f"{stem}_nominal_trajectories.csv", [dict(t=float(tt), arm=cases[j]["arm"], policy=cases[j]["policy"],
                  x=float(y[i,j,0]),s=float(y[i,j,1]),d=float(y[i,j,2])) for i,tt in enumerate(t) for j in selected])


def phase_scan(out, budget, quick):
    rows=[]
    for a in np.linspace(.001,.16,3 if quick else 31):
        for b in np.linspace(0,6,3 if quick else 41):
            budget.check()
            for r in np.linspace(.5,1.8,3 if quick else 41):
                points, continuum = roots("M1",a,b,r)
                kinds=[p["stability"] for p in points]
                rows.append(dict(a=float(a),b=float(b),r=float(r),physical_root_count=len(points),
                                 stable_count=kinds.count("stable"),unstable_count=kinds.count("unstable"),
                                 nonhyperbolic_count=kinds.count("nonhyperbolic"),continuum=continuum,
                                 bistable=kinds==["stable","unstable","stable"],
                                 maximum_residual=max((p["residual"] for p in points),default=0.)))
    csv_write(out/"analytic_phase.csv",rows)
    return dict(points=len(rows),bistable=sum(r["bistable"] for r in rows))


def core_checks(out, budget, quick):
    sys.path.insert(0,str(ROOT/"src"))
    from recoverability.model import Parameters,Schedule,rhs,solve,equilibria,exact_m0
    records=[]
    for model,tau in MODELS[:2] if quick else MODELS:
        p=Parameters(.06,4,1.05,tau=tau)
        schedule=Schedule(.14,.06,1.05)
        state=np.array([.4,.08,.15])
        own=vector_rhs(model,.06,4,1.05,tau,0.,0.,.25,state[None,:])[0]
        err=float(np.max(abs(own-np.asarray(rhs(model,60,state,p)))))
        records.append(dict(check="independent_rhs",model=model,tau=tau,error=err,passed=err<1e-12))
        core=equilibria(model,p)
        ownroots,_=roots(model,.06,4,1.05)
        rooterr=max((abs(q.x-v["x"]) for q,v in zip(core,ownroots)),default=0.)
        records.append(dict(check="independent_roots",model=model,tau=tau,error=rooterr,passed=len(core)==len(ownroots) and rooterr<1e-7))
        for arm,policy,reduce,amp,withdrawal in ARMS[:1] if quick else ARMS:
            budget.check()
            sch=Schedule(.14,.06,1.05,reduce_drive=reduce,repair_amplitude=amp,withdrawal=withdrawal)
            xp=initial_x(model,setting())
            te=np.unique(np.r_[np.linspace(0,120,481),30.,45.])
            reference=solve(model,p,sch,x0=xp,t_eval=te,method="DOP853")
            for dt in (.02,.01,.005) if not quick else (.02,):
                budget.check()
                tested=solve(model,p,sch,x0=xp,t_eval=te,method="RK4",dt=dt)
                delta=float(np.max(abs(tested.y-reference.y)))
                records.append(dict(check="core_rk4_vs_adaptive",model=model,tau=tau,arm=arm,policy=policy,dt=dt,error=delta,passed=delta<1e-6))
                csv_write(out/"numerical_checks.csv",records)
            case=make_cases(model,tau,[setting()])[list(ARMS).index((arm,policy,reduce,amp,withdrawal))]
            tt,yy=integrate(model,[case],120,budget)
            shared=np.searchsorted(tt,te)
            delta=float(np.max(abs(yy[shared,0,:]-reference.y)))
            records.append(dict(check="vectorized_vs_core",model=model,tau=tau,arm=arm,policy=policy,error=delta,passed=delta<1e-6))
        if model=="M0":
            test=solve(model,p,None,x0=.6,t_eval=np.linspace(0,120,481))
            delta=float(np.max(abs(test.x-exact_m0(test.t,.6,.06,1.05))))
            records.append(dict(check="M0_exact",model=model,tau=tau,error=delta,passed=delta<1e-6))
    csv_write(out/"numerical_checks.csv",records)
    if not all(r["passed"] for r in records):
        raise RuntimeError("Stage A numerical checks failed; see numerical_checks.csv")
    return dict(checks=len(records),passed=True)


def treatment_suite(out, budget, quick, parameters):
    allrows=[]
    for model,tau in MODELS:
        cases=make_cases(model,tau,parameters)
        t,y=integrate(model,cases,120,budget)
        rows=outcome_rows(model,cases,t,y,120)
        allrows.extend(rows)
        csv_write(out/"primary_outcomes.csv",allrows)
        save_nominal(out,f"{model}_tau{tau:g}",cases,t,y)
        unresolved=[p for p,row in zip(cases,rows) if row["convergence"]=="unresolved_at_horizon"]
        extended=[]
        for horizon in (240.,480.):
            if not unresolved or quick:
                break
            tt,yy=integrate(model,unresolved,horizon,budget)
            nextrows=outcome_rows(model,unresolved,tt,yy,horizon)
            extended.extend(nextrows)
            unresolved=[p for p,row in zip(unresolved,nextrows) if row["convergence"]=="unresolved_at_horizon"]
        csv_write(out/f"{model}_tau{tau:g}_convergence.csv",extended)
    initial=[]
    for model,tau in MODELS:
        for x0 in (.04,.4,.6,.8):
            cases=make_cases(model,tau,[setting()],x_override=x0)
            t,y=integrate(model,cases,120,budget)
            initial.extend(outcome_rows(model,cases,t,y,120))
    csv_write(out/"initial_state_sensitivity.csv",initial)
    onset=[]
    for model,tau in MODELS:
        for x0 in (.04,.4,.6,.8):
            cases=make_cases(model,tau,[setting()],x_override=x0)
            for case in cases:
                case.update(initial_state_reference="treatment_onset",s0=.14)
            t,y=integrate(model,cases,120,budget)
            onset.extend(outcome_rows(model,cases,t,y,120))
    csv_write(out/"onset_state_sensitivity.csv",onset)
    contrasts=[]
    for model,tau in MODELS:
        for p in parameters:
            subset=[r for r in allrows if r["model"]==model and r["tau"]==tau and r["setting_id"]==p["setting_id"]]
            no=next(r for r in subset if r["arm"]=="none")
            repair=next(r for r in subset if r["arm"]=="repair")
            for policy in ("maintenance","full_withdrawal"):
                drive=next(r for r in subset if r["arm"]=="drive" and r["policy"]==policy)
                both=next(r for r in subset if r["arm"]=="both" and r["policy"]==policy)
                contrasts.append(dict(model=model,tau=tau,setting_id=p["setting_id"],policy=policy,
                    impairment_difference=drive["x_endpoint"]-both["x_endpoint"],
                    target_difference=int(both["primary_target"])-int(drive["primary_target"]),
                    factorial_x=both["x_endpoint"]-drive["x_endpoint"]-repair["x_endpoint"]+no["x_endpoint"],
                    interpretation="synthetic_normalized_design_not_synergy"))
    csv_write(out/"factorial_contrasts.csv",contrasts)
    aggregates=[]
    for model,tau in MODELS:
        for policy in ("maintenance","full_withdrawal"):
            subset=[row for row in contrasts if row["model"]==model and row["tau"]==tau and row["policy"]==policy and row["setting_id"]!="nominal"]
            if subset:
                difference=np.array([row["impairment_difference"] for row in subset])
                count=sum(row["target_difference"]==1 for row in subset)
                aggregates.append(dict(model=model,tau=tau,policy=policy,n_settings=len(subset),
                    positive_target_transition_count=count,design_fraction=count/len(subset),
                    difference_q10=float(np.quantile(difference,.1)),difference_median=float(np.median(difference)),difference_q90=float(np.quantile(difference,.9)),
                    interpretation="arbitrary_uniform_design_not_population_probability; direction_imposed_by_model",
                    ignored_b=model in ("M0","M3")))
    csv_write(out/"design_aggregates.csv",aggregates)
    return dict(primary_cases=len(allrows),initial_sensitivity_cases=len(initial),onset_sensitivity_cases=len(onset),contrasts=len(contrasts))


def alternatives(out,budget,quick):
    structural=[]
    structural_cases=[]
    for d0 in (0.,.15,.30):
        for epsilon in (0.,.003,.01):
            for eta in (0.,.003,.01):
                p=setting(setting_id="nominal")
                p.update(d0=d0,epsilon=epsilon,eta=eta)
                structural_cases.extend(make_cases("M4",10.,[p]))
    t,y=integrate("M4",structural_cases,120,budget)
    structural=outcome_rows("M4",structural_cases,t,y,120)
    csv_write(out/"structural_outcomes.csv",structural)
    unresolved=[p for p,row in zip(structural_cases,structural) if not row["convergence"].startswith("converged")]
    structural_extended=[]
    for horizon in (240.,480.) if not quick else ():
        if not unresolved:
            break
        tt,yy=integrate("M4",unresolved,horizon,budget)
        rows=outcome_rows("M4",unresolved,tt,yy,horizon)
        structural_extended.extend(rows)
        csv_write(out/"structural_convergence.csv",structural_extended)
        unresolved=[p for p,row in zip(unresolved,rows) if not row["convergence"].startswith("converged")]
    costs=[]
    for amplitude in (0.,.2,.4):
        for duration in (5.,15.,30.):
            cases=[p for p in make_cases("M1",10.,[setting(setting_id="nominal")]) if p["arm"]=="both"]
            for p in cases:
                p.update(amplitude=amplitude,pulse_end=30+duration)
                # Full withdrawal remains t=45 even for a planned 30-unit pulse.
            t,y=integrate("M1",cases,120,budget)
            costs.extend(outcome_rows("M1",cases,t,y,120))
    csv_write(out/"amplitude_duration_cost.csv",costs)
    histories=[]
    loops=[]
    for model,tau in MODELS:
        cases=make_cases(model,tau,[setting(setting_id="nominal")])[:2]
        for p,x0 in zip(cases,(.04,.8)):
            p.update(fixed_command=.06,reduce_drive=False,amplitude=0.,x0=x0)
        t,y=integrate(model,cases,480 if not quick else 120,budget)
        for row in outcome_rows(model,cases,t,y,float(t[-1])):
            row.update(experiment="two_histories",final_fixed_a=.06)
            histories.append(row)
        for duration in (60.,240.,960.) if not quick else (60.,):
            p=setting(setting_id="nominal")
            p.update(model=model,tau=tau,arm="ramp",policy="continuous_command",reduce_drive=False,amplitude=0.,withdrawal=None,
                     pulse_end=45.,x0=min(q["x"] for q in roots(model,.02,4.,1.05)[0] if q["stability"]=="stable"))
            tt,yy=integrate(model,[p],2*duration,budget,ramp=duration)
            command=.02+.12*np.where(tt<=duration,tt/duration,2-tt/duration)
            x=yy[:,0,0]
            # Operational crossings are descriptive; they are not fold estimates.
            upward=np.where((tt[1:]<=duration)&(x[:-1]<.4)&(x[1:]>=.4))[0]+1
            downward=np.where((tt[:-1]>=duration)&(x[:-1]>.4)&(x[1:]<=.4))[0]+1
            loops.append(dict(model=model,tau=tau,leg_duration=duration,signed_loop_integral=float(np.trapezoid(x,command)),
                 absolute_enclosed_area=abs(float(np.trapezoid(x,command))),
                 upward_xpoint4_command=float(command[upward[0]]) if len(upward) else None,
                 downward_xpoint4_command=float(command[downward[0]]) if len(downward) else None,
                 initial_above_xpoint4=bool(x[0]>.4),
                 transition_definition="first_directional_x=.4_sample_crossing_not_fold",initial_selection="low_equilibrium_at_a=.02"))
            csv_write(out/f"ramp_{model}_tau{tau:g}_leg{duration:g}.csv",[dict(t=float(a),a_command=float(c),x=float(v),s=float(s)) for a,c,v,s in zip(tt,command,x,yy[:,0,1])])
    csv_write(out/"two_history_endpoints.csv",histories)
    csv_write(out/"ramp_summary.csv",loops)
    return dict(structural_cases=len(structural),cost_cases=len(costs),history_cases=len(histories),ramps=len(loops),
                evolving_structure_root_certification="deferred",M4_convergence_extension_rows=len(structural_extended))


def plots(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    path=out/"primary_outcomes.csv"
    if not path.exists():
        return
    with path.open(encoding="utf-8") as f:
        rows=[r for r in csv.DictReader(f) if r["setting_id"]=="nominal"]
    fig,ax=plt.subplots(figsize=(10,5))
    labels=[f"{r['model']} τ{r['tau']} {r['arm']}\n{r['policy']}" for r in rows]
    ax.bar(np.arange(len(rows)),[float(r["x_endpoint"]) for r in rows])
    ax.axhline(.15,color="black",linestyle="--",label="Synthetic x target .15")
    ax.set_xticks(np.arange(len(rows)),labels,rotation=90,fontsize=6)
    ax.set_ylabel("Dimensionless impairment at t=120")
    ax.set_title("Nominal synthetic model outcomes; no patient interpretation")
    ax.legend();fig.tight_layout()
    fig.savefig(out/"nominal_outcomes.png",dpi=180);fig.savefig(out/"nominal_outcomes.pdf");plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--quick",action="store_true",help="Technical smoke; reduced design, NOT scientific protocol execution")
    mode.add_argument("--full",action="store_true",help="Prespecified A/B/C driver; secondary deferred portions disclosed")
    parser.add_argument("--output-dir",type=Path,required=True,help="New directory only; existing directories rejected")
    parser.add_argument("--stages",choices=("A","B","C","ABC"),default="ABC")
    parser.add_argument("--job-seconds",type=float,default=600.,help="Maximum per-stage budget; cannot exceed 600")
    parser.add_argument("--total-seconds",type=float,default=1800.,help="Overall initial-run budget; cannot exceed 1800")
    args=parser.parse_args(argv)
    if not 0 < args.job_seconds <= 600 or not 0 < args.total_seconds <= 1800:
        parser.error("Budget must be positive and at most 600 seconds/job, 1800 seconds total")
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    parameters=design(args.quick)
    csv_write(out/"parameter_design.csv",parameters)
    csv_write(out/"control_design.csv",[dict(arm=a,policy=p,reduce_drive=d,repair_amplitude=r,withdrawal=w,onset=30,pulse_end=45) for a,p,d,r,w in ARMS])
    inputs=[ROOT/"protocol/study_design_v1.md",Path(__file__),ROOT/"src/recoverability/model.py"]
    if (ROOT/"protocol/clarifications_v1_1.md").exists():
        inputs.append(ROOT/"protocol/clarifications_v1_1.md")
    metadata=dict(evidence_class="technical_smoke" if args.quick else "synthetic_simulation",
                  protocol_complete=False,stages_requested=args.stages,seed=SEED,dimensionless=True,
                  rtol=1e-9,atol=1e-11,max_step=MAX_STEP,solver="DOP853",sample_step=SAMPLE_STEP,
                  code_protocol_hashes={str(p.relative_to(ROOT)):sha(p) if p.exists() else "missing" for p in inputs},
                  python=platform.python_version(),packages={p:importlib.metadata.version(p) for p in ("numpy","scipy","matplotlib")},
                  owner="literature_experiment_driver",deferred=["Stage D observations/nuisance/identifiability/fitting",
                  "M4 evolving-structure isolated/manifold equilibrium certification",
                  "Global control-set reachability proof; pulse-family results only",
                  "Continuous persistence extrema certification; sampled maxima reported"],jobs=[],failures=[])
    (out/"manifest.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    jobs={"A":lambda b:dict(core=core_checks(out,b,args.quick),phase=phase_scan(out,b,args.quick)),
          "B":lambda b:treatment_suite(out,b,args.quick,parameters),
          "C":lambda b:alternatives(out,b,args.quick)}
    for stage in args.stages:
        before=time.monotonic()
        job=dict(stage=stage,owner="literature_experiment_driver",inputs=["protocol", "core model", "parameter/control tables"],
                 dependency="Stage A verified core for B/C",acceptance="Numerical checks; bounds; no clipping; explicit unresolved states",status="running")
        try:
            if stage in ("B","C") and args.stages=="ABC" and any(j["status"]!="generated_unreviewed" for j in metadata["jobs"]):
                raise RuntimeError("Earlier stage failure prevents dependent stage")
            remaining=args.total_seconds-(time.monotonic()-start)
            if remaining<=0:
                raise TimeoutError("Overall run budget exhausted")
            result=jobs[stage](Budget(min(args.job_seconds,remaining)))
            job.update(status="generated_unreviewed",result=result)
        except Exception as exc:
            job.update(status="failed_or_incomplete",error=str(exc))
            metadata["failures"].append(dict(stage=stage,type=type(exc).__name__,message=str(exc)))
            (out/f"error_stage_{stage}.txt").write_text(traceback.format_exc(),encoding="utf-8")
        job["elapsed_seconds"]=time.monotonic()-before
        metadata["jobs"].append(job)
        (out/"summary.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
        if job["status"]!="generated_unreviewed":
            break
    try:
        plots(out)
    except Exception as exc:
        metadata["failures"].append(dict(stage="plots",message=str(exc)))
    metadata["elapsed_seconds"]=time.monotonic()-start
    metadata["status"]="generated_partial_protocol_unreviewed" if not metadata["failures"] else "failed_or_incomplete"
    (out/"summary.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    (out/"output_hashes.json").write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="output_hashes.json"},indent=2),encoding="utf-8")
    print(json.dumps(dict(status=metadata["status"],output_dir=str(out),failures=metadata["failures"])))
    return 1 if metadata["failures"] else 0


if __name__=="__main__":
    raise SystemExit(main())
