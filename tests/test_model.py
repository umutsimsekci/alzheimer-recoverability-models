"""Mathematical and event-handling checks, not biological validation."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from recoverability.model import Parameters, Schedule, equilibria, exact_m0, rhs, solve


def test_nominal_roots_and_stability():
    points = equilibria("M1", Parameters(.06, 4, 1.05))
    np.testing.assert_allclose([p.x for p in points], [.070871215252208, .4, .529128784747792], atol=1e-12)
    assert [p.stability for p in points] == ["stable", "unstable", "stable"]
    assert max(p.residual for p in points) < 1e-12


def test_m2_has_one_physical_stable_root_and_m3_matches_m0():
    p = Parameters(.06, 4, 1.05)
    roots = equilibria("M2", p)
    assert len(roots) == 1 and roots[0].stability == "stable"
    assert equilibria("M3", p) == equilibria("M0", p)


def test_fold_and_cusp_are_not_hyperbolic_stable_points():
    x, b = .21, 4.
    roots = equilibria("M1", Parameters(b*x*x*(1-2*x), b, 2*b*x*(1-x)**2))
    assert len(roots) == 2
    assert sum(p.stability == "nonhyperbolic" for p in roots) == 1
    cusp = equilibria("M1", Parameters(b/27, b, 8*b/27))
    assert len(cusp) == 1 and cusp[0].stability == "nonhyperbolic"
    assert cusp[0].x == pytest.approx(1/3, abs=1e-10)


def test_boundary_roots_and_continuum():
    zero = Parameters(0, 0, 0)
    for model in ("M0", "M1", "M2", "M3", "M4"):
        assert equilibria(model, zero).continuum
    roots = equilibria("M1", Parameters(0, 4, 1))
    assert roots[0].x == 0 and roots[0].stability == "stable"
    assert any(p.stability == "nonhyperbolic" for p in roots)
    assert equilibria("M1", Parameters(.06, 4, 0))[-1].x == 1


@pytest.mark.parametrize("method", ["DOP853", "RK4"])
def test_m0_against_exact_solution_and_zero_drift(method):
    times = np.linspace(0, 20, 201)
    p = Parameters(.06, 0, 1.05)
    result = solve("M0", p, x0=.5, t_eval=times, method=method)
    np.testing.assert_allclose(result.x, exact_m0(times, .5, .06, 1.05), atol=2e-9, rtol=0)
    constant = solve("M0", Parameters(0, 0, 0), x0=.6, d0=.3, t_eval=times, method=method)
    np.testing.assert_array_equal(constant.x, np.full(len(times), .6))
    np.testing.assert_array_equal(constant.d, np.full(len(times), .3))


def test_exact_segmentation_with_off_grid_events_and_withdrawal():
    times = np.array([0, .2, .9, 1.8, 2.5, 3., 4.])
    schedule = Schedule(.14, .06, 1.05, onset=.37, pulse_end=1.43, withdrawal=2.17)
    expected = []
    for target in times:
        value, start = .6, 0.
        boundaries = sorted({0., target, *[e for e in schedule.events() if 0 < e < target]})
        for end in boundaries[1:]:
            a, r = schedule.coefficients((start+end)/2)
            value = float(exact_m0(end, value, a, r, t0=start))
            start = end
        expected.append(value)
    adaptive = solve("M0", Parameters(.06, 0, 1.05), schedule, .6, t_eval=times)
    explicit = solve("M0", Parameters(.06, 0, 1.05), schedule, .6, t_eval=times, method="RK4", dt=.02)
    assert adaptive.metadata["events"] == [.37, 1.43, 2.17]
    assert adaptive.metadata["segments"] == [[0., .37], [.37, 1.43], [1.43, 2.17], [2.17, 4.]]
    np.testing.assert_allclose(adaptive.x, expected, atol=2e-9, rtol=0)
    np.testing.assert_allclose(explicit.x, expected, atol=2e-9, rtol=0)
    np.testing.assert_array_equal(adaptive.t, times)


def test_schedule_right_continuity_and_early_withdrawal():
    schedule = Schedule(.14, .06, 1.05, onset=1., pulse_end=3., withdrawal=2.)
    assert schedule.coefficients(.999) == (.14, 1.05)
    assert schedule.coefficients(1.) == pytest.approx((.06, 1.45))
    assert schedule.coefficients(2.) == (.14, 1.05)
    assert schedule.coefficients(3.) == (.14, 1.05)


def test_nominal_maintenance_vs_full_withdrawal_and_rk4_refinement():
    p = Parameters(.06, 4, 1.05)
    times = np.linspace(0, 120, 241)
    drive = Schedule(.14, .06, 1.05, repair_amplitude=0)
    combined = Schedule(.14, .06, 1.05)
    lower = solve("M1", p, drive, .6, t_eval=times)
    maintained = solve("M1", p, combined, .6, t_eval=times)
    withdrawn = solve("M1", p, Schedule(.14, .06, 1.05, withdrawal=45), .6, t_eval=times)
    assert lower.x[-1] == pytest.approx(.529128784747792, abs=1e-8)
    assert maintained.x[-1] == pytest.approx(.070871215252208, abs=1e-8)
    assert withdrawn.x[-1] == pytest.approx(.602781900257604, abs=1e-8)
    assert np.all(maintained.x <= lower.x + 1e-8)
    for dt in (.02, .01, .005):
        independent = solve("M1", p, combined, .6, t_eval=times, method="RK4", dt=dt)
        np.testing.assert_allclose(independent.x, maintained.x, atol=1e-7, rtol=0)


def test_m3_slow_state_analytic_response_and_steady_linear_limit():
    p = Parameters(.06, 4, 1.05, tau=10)
    schedule = Schedule(.14, .06, 1.05, onset=1., pulse_end=1., repair_amplitude=0)
    times = np.linspace(0, 40, 161)
    result = solve("M3", p, schedule, .4, t_eval=times)
    expected = np.where(times <= 1, .14, .06+.08*np.exp(-(times-1)/10))
    np.testing.assert_allclose(result.s, expected, atol=2e-10, rtol=0)
    frozen = solve("M3", p, x0=.4, t_eval=times)
    np.testing.assert_allclose(frozen.x, exact_m0(times, .4, .06, 1.05), atol=2e-9, rtol=0)


def test_m4_bounds_frozen_structure_and_regeneration():
    p = Parameters(.06, 4, 1.05, epsilon=.01, eta=.003)
    result = solve("M4", p, x0=.6, d0=.3, t_eval=np.linspace(0, 120, 241))
    assert result.x.min() >= 0 and result.x.max() <= 1
    assert result.d.min() >= 0 and result.d.max() <= 1
    fixed = Parameters(.06, 4, 1.05)
    roots = equilibria("M4", fixed, d=.3)
    transformed = equilibria("M1", Parameters(.06, 4, 1.05*.7))
    np.testing.assert_allclose([p.x for p in roots], [p.x for p in transformed])
    no_drive = Parameters(0, 0, 1, epsilon=0, eta=.1)
    times = np.linspace(0, 20, 81)
    restored = solve("M4", no_drive, x0=0, d0=.5, t_eval=times)
    np.testing.assert_allclose(restored.d, .5*np.exp(-.1*times), atol=2e-9, rtol=0)
    complete_loss = solve("M4", fixed, x0=1, d0=1, t_eval=[0., 1.])
    np.testing.assert_array_equal(complete_loss.x, [1., 1.])
    np.testing.assert_array_equal(complete_loss.d, [1., 1.])


def test_m4_boundary_derivatives_and_threshold():
    p = Parameters(.06, 4, 1.05, epsilon=.01, eta=.003)
    assert rhs("M4", 0, [0, .06, 0], p)[0] >= 0
    assert rhs("M4", 0, [1, .06, 1], p)[0] <= 0
    assert rhs("M4", 0, [.6, .06, 0], p)[2] >= 0
    assert rhs("M4", 0, [.6, .06, 1], p)[2] <= 0
    assert rhs("M4", 0, [p.xc, .06, 0], p)[2] == 0


@pytest.mark.parametrize("kwargs", [{"a": -1}, {"r": np.nan}, {"tau": 0}, {"xc": 1.1}])
def test_invalid_parameters(kwargs):
    fields = dict(a=.06, b=4., r=1.05)
    fields.update(kwargs)
    with pytest.raises(ValueError):
        Parameters(**fields)


@pytest.mark.parametrize("kwargs", [{"pulse_end": 20}, {"withdrawal": 20}, {"repair_amplitude": -1}, {"a_post": .2}])
def test_invalid_schedules(kwargs):
    fields = dict(a_pre=.14, a_post=.06, r_base=1.05)
    fields.update(kwargs)
    with pytest.raises(ValueError):
        Schedule(**fields)


@pytest.mark.parametrize("times", [[], [0, 0], [2, 1], [-1, 0], [0, np.inf], [[0, 1]]])
def test_invalid_timegrids(times):
    with pytest.raises(ValueError):
        solve("M1", Parameters(.06, 4, 1.05), t_eval=times)


def test_single_time_initial_state_and_invalid_integrator():
    p = Parameters(.06, 4, 1.05)
    single = solve("M1", p, x0=.6, t_eval=[2.])
    np.testing.assert_array_equal(single.y, [[.6, .06, 0.]])
    with pytest.raises(ValueError):
        solve("M1", p, method="Euler")
    with pytest.raises(ValueError):
        solve("M1", p, dt=0)
    with pytest.raises(ValueError):
        solve("M9", p)


def test_coarse_rk4_violation_raises_without_clipping():
    with pytest.raises(ArithmeticError, match="bounds"):
        solve("M0", Parameters(0, 0, 10), x0=.5, t_eval=[0, 1], method="RK4", dt=1)
