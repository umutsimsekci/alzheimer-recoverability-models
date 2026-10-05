"""Independent analytic and per-case validation of the stacked study driver."""
from pathlib import Path
import importlib.util
from types import SimpleNamespace
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("study_driver_validation", ROOT / "scripts/run_study_v1.py")
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
sys.path.insert(0, str(ROOT / "src"))
from recoverability.model import Parameters, Schedule, solve


def exact_scalar(model, elapsed, x0, a, b, r):
    """Closed-form M0 exponential and M2 Riccati, independent of core code."""
    if model == "M0" or b == 0:
        rate = a+r
        return x0 if rate == 0 else a/rate+(x0-a/rate)*np.exp(-rate*elapsed)
    linear = b-a-r
    gap = np.sqrt(linear*linear+4*b*a)
    upper, lower = (linear+gap)/(2*b), (linear-gap)/(2*b)
    ratio = (x0-upper)/(x0-lower)*np.exp(-gap*elapsed)
    return (upper-ratio*lower)/(1-ratio)


@pytest.mark.parametrize("model", ["M0", "M2"])
def test_offgrid_eventwise_analytic_reference(model):
    case = driver.make_cases(model, 10, [driver.setting()])[5]
    case.update(pulse_end=42.37, withdrawal=41.13, x0=.6)
    t, y = driver.integrate(model, [case], 120, driver.Budget(30))
    expected=[]
    for target in t:
        cuts=sorted({0.,float(target),*[e for e in (30.,41.13,42.37) if 0<e<target]})
        value=.6
        for left,right in zip(cuts[:-1],cuts[1:]):
            middle=(left+right)/2
            a=.06 if 30<=middle<41.13 else .14
            r=1.45 if 30<=middle<min(42.37,41.13) else 1.05
            value=exact_scalar(model,right-left,value,a,4.,r)
        expected.append(value)
    np.testing.assert_allclose(y[:,0,0],expected,atol=1e-7,rtol=0)


def test_m3_analytic_drive_and_m4_analytic_regeneration():
    case=driver.make_cases("M3",50,[driver.setting()])[2]
    t,y=driver.integrate("M3",[case],120,driver.Budget(30))
    expected=np.where(t<30,.14,.06+.08*np.exp(-(t-30)/50))
    np.testing.assert_allclose(y[:,0,1],expected,atol=1e-9,rtol=0)
    case=driver.make_cases("M4",10,[driver.setting()])[0]
    case.update(fixed_command=0.,x0=0.,d0=.5,epsilon=0.,eta=.1)
    t,y=driver.integrate("M4",[case],120,driver.Budget(30))
    np.testing.assert_allclose(y[:,0,2],.5*np.exp(-.1*t),atol=1e-9,rtol=0)


def test_controls_cost_and_policy_deduplication():
    cases=driver.make_cases("M1",10,[driver.setting()])
    assert len(cases)==6
    assert [(p["arm"],p["policy"]) for p in cases[:2]]==[("none","identical_policies"),("repair","identical_policies")]
    t,y=driver.integrate("M1",cases,120,driver.Budget(30))
    rows=driver.outcome_rows("M1",cases,t,y,120)
    np.testing.assert_allclose([r["cost_total"] for r in rows],[0,15,90,15,105,30],atol=1e-12)
    assert rows[4]["x_endpoint"]==pytest.approx(.070871215252208,abs=1e-8)
    assert rows[5]["x_endpoint"]==pytest.approx(.602781900257604,abs=1e-8)
    assert rows[4]["primary_target"] and not rows[5]["primary_target"]
    cases[5]["pulse_end"]=60.
    t,y=driver.integrate("M1",[cases[5]],120,driver.Budget(30))
    assert driver.outcome_rows("M1",[cases[5]],t,y,120)[0]["cost_total"]==30


def test_convergence_checks_history_drift_and_hidden_state():
    case=driver.make_cases("M0",10,[driver.setting()])[2]
    t=np.linspace(0,120,481)
    equilibrium=.06/(.06+1.05)
    y=np.zeros((len(t),1,3));y[:,:,0]=equilibrium;y[:,:,1]=.14
    y[(t>=110)&(t<120),0,0]+= .1
    assert driver.outcome_rows("M0",[case],t,y,120)[0]["convergence"]=="unresolved_at_horizon"
    y[:,:,0]=equilibrium;y[:,:,1]=.2
    assert driver.outcome_rows("M3",[case],t,y,120)[0]["convergence"]=="unresolved_at_horizon"
    case.update(epsilon=.003,eta=0.,d0=1.)
    y[:,:,0]=1.;y[:,:,2]=1.
    row=driver.outcome_rows("M4",[case],t,y,120)[0]
    assert row["convergence"]=="full_drift_small_equilibrium_certification_deferred"
    assert row["equilibrium_flag_scope"]=="frozen_d_only"


@pytest.mark.parametrize("state_column,bad_value",[(0,1.1),(2,-.1),(1,-.1)])
def test_internal_accepted_steps_are_checked(monkeypatch,state_column,bad_value):
    case=driver.make_cases("M1",10,[driver.setting()])[0]
    def counterfeit(fun,span,state,**kwargs):
        hidden=np.asarray(state).copy();hidden[state_column]=bad_value
        return SimpleNamespace(success=True,y=np.column_stack([hidden,state]),
            sol=lambda times:np.repeat(np.asarray(state)[:,None],len(times),axis=1))
    monkeypatch.setattr(driver,"solve_ivp",counterfeit)
    with pytest.raises(RuntimeError,match="bounds|nonnegative|Nonfinite|negative"):
        driver.integrate("M1",[case],120,driver.Budget(30))


def test_full_m1_stack_against_individual_cases_at_120_and_480():
    parameters=driver.design(False)
    cases=driver.make_cases("M1",10,parameters)
    t,y=driver.integrate("M1",cases,480,driver.Budget(90))
    selected=[0,1,2,8,128,256]  # nominal + sample indices 0,1,7,127,255
    for setting_index in selected:
        for arm_index in range(6):
            index=setting_index*6+arm_index
            case=cases[index]
            schedule=Schedule(case["a"]+.08,case["a"],case["r"],
                pulse_end=case["pulse_end"],withdrawal=case["withdrawal"],
                reduce_drive=case["reduce_drive"],repair_amplitude=case["amplitude"])
            individual=solve("M1",Parameters(case["a"],case["b"],case["r"]),
                schedule,case["x0"],t_eval=t)
            np.testing.assert_allclose(y[:,index],individual.y,atol=1e-6,rtol=0)
            for endpoint in (120.,480.):
                position=np.flatnonzero(t==endpoint)[0]
                assert np.max(abs(y[position,index]-individual.y[position]))<1e-6


@pytest.mark.parametrize("model,tau",driver.MODELS)
def test_nominal_and_stress_batch_vs_individual(model,tau):
    fold=.1027147093197253
    parameters=[driver.setting(),driver.setting(a=fold-1e-5),driver.setting(a=fold+1e-5)]
    cases=driver.make_cases(model,tau,parameters)
    t,y=driver.integrate(model,cases,480,driver.Budget(45))
    for index,case in enumerate(cases):
        schedule=Schedule(case["a"]+.08,case["a"],case["r"],withdrawal=case["withdrawal"],
            reduce_drive=case["reduce_drive"],repair_amplitude=case["amplitude"])
        independent=solve(model,Parameters(case["a"],case["b"],case["r"],tau=tau),
            schedule,case["x0"],t_eval=t)
        np.testing.assert_allclose(y[:,index],independent.y,atol=1e-6,rtol=0)


def test_two_history_analytic_m0_and_distinct_m1_equilibria():
    for model in ("M0","M1","M2","M3"):
        cases=driver.make_cases(model,50,[driver.setting()])[:2]
        for case,x0 in zip(cases,(.04,.8)):
            case.update(fixed_command=.06,reduce_drive=False,amplitude=0.,x0=x0)
        t,y=driver.integrate(model,cases,480,driver.Budget(30))
        if model=="M1":
            np.testing.assert_allclose(y[-1,:,0],[.070871215252208,.529128784747792],atol=1e-8)
        else:
            assert abs(y[-1,0,0]-y[-1,1,0])<1e-8
        if model=="M0":
            for index,x0 in enumerate((.04,.8)):
                np.testing.assert_allclose(y[:,index,0],exact_scalar("M0",t,x0,.06,0,1.05),atol=1e-8)


def test_m3_continuous_ramp_drive_against_closed_form():
    duration,tau=60.,50.
    case=driver.make_cases("M3",tau,[driver.setting()])[0]
    case.update(x0=.02/(.02+1.05))
    t,y=driver.integrate("M3",[case],2*duration,driver.Budget(30),ramp=duration)
    slope=.12/duration
    first=.02+slope*(np.minimum(t,duration)-tau)+slope*tau*np.exp(-np.minimum(t,duration)/tau)
    at_turn=.02+slope*(duration-tau)+slope*tau*np.exp(-duration/tau)
    elapsed=np.maximum(t-duration,0)
    second=.14-slope*(elapsed-tau)+(at_turn-.14-slope*tau)*np.exp(-elapsed/tau)
    expected=np.where(t<=duration,first,second)
    np.testing.assert_allclose(y[:,0,1],expected,atol=1e-9,rtol=0)


def test_m4_frozen_and_evolving_batches_against_individual_cases():
    cases=[]
    for d0,epsilon,eta in ((.15,0.,0.),(.3,.003,0.),(.3,.003,.01)):
        parameters=driver.setting()
        parameters.update(d0=d0,epsilon=epsilon,eta=eta)
        cases.extend(driver.make_cases("M4",10,[parameters]))
    t,y=driver.integrate("M4",cases,480,driver.Budget(45))
    for index,case in enumerate(cases):
        schedule=Schedule(case["a"]+.08,case["a"],case["r"],withdrawal=case["withdrawal"],
            reduce_drive=case["reduce_drive"],repair_amplitude=case["amplitude"])
        independent=solve("M4",Parameters(case["a"],case["b"],case["r"],epsilon=case["epsilon"],eta=case["eta"]),
            schedule,case["x0"],d0=case["d0"],t_eval=t)
        np.testing.assert_allclose(y[:,index],independent.y,atol=1e-6,rtol=0)


def test_failure_is_recorded_and_existing_output_not_overwritten(monkeypatch,tmp_path):
    import json
    output=tmp_path/"deliberate_failure"
    def forced_failure(*args):
        raise ArithmeticError("deliberate validation failure")
    monkeypatch.setattr(driver,"core_checks",forced_failure)
    monkeypatch.setattr(driver,"plots",lambda out:None)
    assert driver.main(["--quick","--stages","A","--output-dir",str(output)])==1
    result=json.loads((output/"summary.json").read_text())
    assert result["status"]=="failed_or_incomplete"
    assert result["jobs"][0]["status"]=="failed_or_incomplete"
    assert result["failures"][0]["type"]=="ArithmeticError"
    assert "deliberate validation failure" in (output/"error_stage_A.txt").read_text()
    assert (output/"output_hashes.json").exists()
    with pytest.raises(FileExistsError):
        driver.main(["--quick","--stages","A","--output-dir",str(output)])
