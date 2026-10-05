from dataclasses import replace
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/run_observation_audit.py"
spec = importlib.util.spec_from_file_location("observation_audit", SCRIPT)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_event_times_and_half_open_nuisance_window():
    times = audit.observation_times()
    assert len(times) == 62
    assert {0.,30.,45.,60.,120.}.issubset(set(times))
    assert np.all(np.diff(times)>0)
    np.testing.assert_array_equal(audit.nuisance(np.array([29.999,30.,59.999,60.])),[0.,.12,.12,0.])


def test_observation_changes_only_readout_and_shared_nuisance_is_explicit():
    times = np.array([0.,30.,45.,60.])
    x, d = np.full(4,.4), np.full(4,.15)
    before_x,before_d = x.copy(),d.copy()
    clean,offset,channels = audit.measurement_values(times,x,d,"nuisance_y1_only")
    assert channels==2
    np.testing.assert_array_equal(offset[:,1],0.)
    shared_clean,shared,_ = audit.measurement_values(times,x,d,"nuisance_shared_raw")
    np.testing.assert_array_equal(shared[:,0],shared[:,1])
    np.testing.assert_array_equal(clean,shared_clean)
    np.testing.assert_array_equal(x,before_x); np.testing.assert_array_equal(d,before_d)


def test_instantaneous_rank_distinguishes_redundancy_structure_and_nuisance_loadings():
    rows = {row["design"]:row for row in audit.observation_jacobians()}
    assert rows["anchored_x"]["instantaneous_rank"]==1
    assert rows["redundant_tasks"]["instantaneous_rank"]==1
    assert rows["anchored_x_d"]["instantaneous_rank"]==2
    assert rows["nuisance_shared_raw_unknown_nuisance"]["instantaneous_rank"]==2
    assert rows["nuisance_shared_loading_scaled_algebraic_alternative"]["instantaneous_rank"]==1


def test_design_fixed_seed_and_indices_are_prespecified():
    first,second = audit.parameter_design(),audit.parameter_design()
    assert first.shape==(256,3)
    np.testing.assert_array_equal(first,second)
    assert np.all(first>=[.03,2.,.8]) and np.all(first<=[.09,6.,1.3])
    assert audit.SELECTED_INDICES==tuple(range(8))


def test_logparameter_sensitivity_matches_independent_central_difference():
    p = audit.Parameters(.06,4.,1.05)
    times = audit.observation_times()
    x0 = audit.initial_state(p)
    schedule = audit.make_schedule(p,"combination")
    matrix,error = audit.scalar_sensitivities(p,schedule,x0,times)
    assert error<1e-7
    # Core trajectory absolute errors are amplified by 1/h; 1e-3 avoids
    # subtractive numerical noise while retaining second-order truncation.
    h=1e-3
    for column,name in enumerate(("a","b","r")):
        plus=replace(p,**{name:getattr(p,name)*np.exp(h)})
        minus=replace(p,**{name:getattr(p,name)*np.exp(-h)})
        y_plus=1.-audit.solve("M1",plus,audit.make_schedule(plus),x0=x0,t_eval=times).x
        y_minus=1.-audit.solve("M1",minus,audit.make_schedule(minus),x0=x0,t_eval=times).x
        np.testing.assert_allclose(matrix[:,column],(y_plus-y_minus)/(2*h),rtol=3e-4,atol=3e-5)


def test_steady_trajectory_has_rank_one_not_three():
    p = audit.Parameters(.06,4.,1.05)
    matrix,_ = audit.scalar_sensitivities(p,audit.make_schedule(p,"baseline"),audit.initial_state(p),audit.observation_times())
    assert audit.rank_record(matrix)["rank"]==1


def test_exact_confounding_and_known_additive_pulse_qualification():
    algebra,_,_ = audit.algebraic_examples(audit.observation_times())
    confounding=algebra["frozen_m4_confounding"]
    assert confounding["max_anchored_x_difference_equivalent_pair"]<1e-7
    assert confounding["max_x_difference_known_r_changed_d"]>1e-3
    assert confounding["max_x_difference_known_additive_pulse_pair"]>1e-3
    assert confounding["rank"]==1
    assert algebra["steady_state_scaling"]["max_root_difference"]<1e-10
    assert algebra["steady_state_scaling"]["max_clock_time_trajectory_difference"]>1e-3
    assert algebra["task_only_false_improvement"]["latent_max_drift_from_initial"]<1e-7
    assert algebra["task_only_false_improvement"]["d_max_drift_from_initial"]==0.


def test_invalid_measurement_design_is_rejected():
    with pytest.raises(ValueError):
        audit.measurement_values(np.array([0.]),np.array([.4]),np.array([.15]),"unknown")
