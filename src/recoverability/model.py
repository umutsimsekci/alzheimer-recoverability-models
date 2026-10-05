"""Deterministic synthetic model core. Time/rates are uncalibrated model units.

All trajectories use [x, s, d]: impairment, slow drive, structural loss.
Only M3 evolves s; only M4 evolves or uses d. No state is clipped.
"""
from dataclasses import dataclass, field
from typing import Iterator

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MODELS = frozenset({"M0", "M1", "M2", "M3", "M4"})
RTOL, ATOL = 1e-9, 1e-11


def _finite(value, name):
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")


def _nonnegative(value, name):
    _finite(value, name)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")


def _model(model_id):
    if model_id not in MODELS:
        raise ValueError(f"Unknown model {model_id!r}; expected M0--M4")


@dataclass(frozen=True)
class Parameters:
    a: float
    b: float
    r: float
    tau: float = 10.0
    epsilon: float = 0.0
    eta: float = 0.0
    xc: float = 0.25

    def __post_init__(self):
        for name in ("a", "b", "r", "epsilon", "eta"):
            _nonnegative(getattr(self, name), name)
        _finite(self.tau, "tau")
        if self.tau <= 0:
            raise ValueError("tau must be positive")
        _finite(self.xc, "xc")
        if not 0 <= self.xc <= 1:
            raise ValueError("xc must lie in [0,1]")


@dataclass(frozen=True)
class Schedule:
    a_pre: float
    a_post: float
    r_base: float
    onset: float = 30.0
    pulse_end: float = 45.0
    withdrawal: float | None = None
    reduce_drive: bool = True
    repair_amplitude: float = 0.4

    def __post_init__(self):
        for name in ("a_pre", "a_post", "r_base", "onset", "pulse_end", "repair_amplitude"):
            _nonnegative(getattr(self, name), name)
        if self.pulse_end < self.onset:
            raise ValueError("pulse_end must not precede onset")
        if self.reduce_drive and self.a_post > self.a_pre:
            raise ValueError("drive reduction requires a_post <= a_pre")
        if self.withdrawal is not None:
            _nonnegative(self.withdrawal, "withdrawal")
            if self.withdrawal < self.onset:
                raise ValueError("withdrawal must not precede onset")

    def coefficients(self, t: float) -> tuple[float, float]:
        """Right-continuous commanded a and effective repair r."""
        _finite(t, "t")
        active = t >= self.onset and (self.withdrawal is None or t < self.withdrawal)
        a = self.a_post if active and self.reduce_drive else self.a_pre
        r = self.r_base + (self.repair_amplitude if active and t < self.pulse_end else 0.0)
        return a, r

    def events(self) -> tuple[float, ...]:
        values = [self.onset, self.pulse_end]
        if self.withdrawal is not None:
            values.append(self.withdrawal)
        return tuple(sorted(set(values)))


@dataclass(frozen=True)
class Equilibrium:
    x: float
    derivative: float
    stability: str
    residual: float


@dataclass(frozen=True)
class EquilibriumSet:
    points: tuple[Equilibrium, ...]
    continuum: bool = False

    def __iter__(self) -> Iterator[Equilibrium]:
        return iter(self.points)

    def __len__(self):
        return len(self.points)

    def __getitem__(self, item):
        return self.points[item]


@dataclass
class Solution:
    t: np.ndarray
    y: np.ndarray
    model_id: str
    method: str
    metadata: dict = field(default_factory=dict)

    @property
    def x(self):
        return self.y[:, 0]

    @property
    def s(self):
        return self.y[:, 1]

    @property
    def d(self):
        return self.y[:, 2]


def _rates(model_id, y, p, a, r):
    x, s, d = y
    if model_id in ("M0", "M3"):
        drive = s if model_id == "M3" else a
        dx = drive * (1 - x) - r * x
    elif model_id == "M2":
        dx = a * (1 - x) + p.b * x * (1 - x) - r * x
    else:
        repair = r * (1 - d) if model_id == "M4" else r
        dx = a * (1 - x) + p.b * x * x * (1 - x) - repair * x
    ds = (a - s) / p.tau if model_id == "M3" else 0.0
    dd = p.epsilon * max(x - p.xc, 0.0) * (1 - d) - p.eta * d if model_id == "M4" else 0.0
    return np.array([dx, ds, dd], dtype=float)


def rhs(model_id, t, y, params, schedule=None):
    """Public derivative, including right-continuous schedules at events."""
    _model(model_id)
    y = np.asarray(y, dtype=float)
    if y.shape != (3,) or not np.isfinite(y).all():
        raise ValueError("y must be a finite vector [x,s,d]")
    _finite(t, "t")
    a, r = schedule.coefficients(t) if schedule is not None else (params.a, params.r)
    return _rates(model_id, y, params, a, r)


def equilibria(model_id, params, d=0.0, residual_tol=1e-10, stability_tol=1e-7):
    """Physical fixed-input x equilibria. M4 freezes d; M3 sets s=a.

    Derivative-zero points are labelled nonhyperbolic, never silently stable.
    A continuum means every x in [0,1] is an equilibrium, not an empty set.
    """
    _model(model_id)
    _finite(d, "d")
    if not 0 <= d <= 1:
        raise ValueError("d must lie in [0,1]")
    for value, name in ((residual_tol, "residual_tol"), (stability_tol, "stability_tol")):
        _finite(value, name)
        if value <= 0:
            raise ValueError(f"{name} must be positive")
    a, b = params.a, params.b
    r = params.r * (1 - d) if model_id == "M4" else params.r
    if model_id in ("M0", "M3"):
        coefficients = np.array([-(a + r), a])
    elif model_id == "M2":
        coefficients = np.array([-b, b - a - r, a])
    else:
        coefficients = np.array([-b, b, -(a + r), a])
    coefficients = np.trim_zeros(coefficients, "f")
    if not coefficients.size:
        return EquilibriumSet((), continuum=True)
    poly = np.poly1d(coefficients)
    deriv = np.polyder(poly)
    # At the cusp the derivative has a repeated real root. Roundoff can give
    # its quadratic roots imaginary parts ~1e-9; certify the resulting
    # tangent with the original drift residual instead of discarding it.
    critical = [float(z.real) for z in np.roots(deriv)
                if abs(z.imag) < 1e-7 and 0 < z.real < 1] if deriv.order > 0 else []
    cuts = sorted(set([0.0, *critical, 1.0]))
    roots = [x for x in cuts if abs(poly(x)) <= residual_tol]
    for left, right in zip(cuts[:-1], cuts[1:]):
        if abs(poly(left)) > residual_tol and abs(poly(right)) > residual_tol and poly(left) * poly(right) < 0:
            roots.append(brentq(poly, left, right, xtol=5e-15, rtol=1e-14))
    unique = []
    for x in sorted(roots):
        if not unique or abs(x - unique[-1]) > 1e-8:
            unique.append(x)
    points = []
    for x in unique:
        derivative = float(deriv(x))
        residual = float(abs(poly(x)))
        if residual > residual_tol:
            raise ArithmeticError("Equilibrium residual exceeds tolerance")
        stability = "stable" if derivative < -stability_tol else "unstable" if derivative > stability_tol else "nonhyperbolic"
        points.append(Equilibrium(float(x), derivative, stability, residual))
    return EquilibriumSet(tuple(points))


def exact_m0(t, x0, a, r, t0=0.0):
    """Exact constant-coefficient M0 trajectory, including zero drift."""
    for value, name in ((a, "a"), (r, "r")):
        _nonnegative(value, name)
    _finite(x0, "x0")
    _finite(t0, "t0")
    if not 0 <= x0 <= 1:
        raise ValueError("x0 must lie in [0,1]")
    t = np.asarray(t, dtype=float)
    if not np.isfinite(t).all() or np.any(t < t0):
        raise ValueError("t must be finite and >= t0")
    rate = a + r
    if rate == 0:
        return np.full_like(t, x0)
    equilibrium = a / rate
    return equilibrium + (x0 - equilibrium) * np.exp(-rate * (t - t0))


def _check_bounds(y):
    if not np.isfinite(y).all():
        raise ArithmeticError("Nonfinite integration state")
    if np.any(y[..., [0, 2]] < -1e-9) or np.any(y[..., [0, 2]] > 1 + 1e-9) or np.any(y[..., 1] < -1e-9):
        raise ArithmeticError("Integration violates physical bounds; no clipping applied")


def _rk4_step(y, h, derivative):
    k1 = derivative(y)
    k2 = derivative(y + .5 * h * k1)
    k3 = derivative(y + .5 * h * k2)
    k4 = derivative(y + h * k3)
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def solve(model_id, params, schedule=None, x0=.04, d0=0.0, s0=None,
          t_eval=None, method="DOP853", dt=.02):
    """Integrate piecewise-constant controls with exact event segmentation.

    Initial state is defined at t_eval[0]. Output times equal t_eval exactly.
    Schedule commands override params.a/r; other coefficients use params.
    """
    _model(model_id)
    if method not in ("DOP853", "RK4"):
        raise ValueError("method must be 'DOP853' or 'RK4'")
    _finite(dt, "dt")
    if dt <= 0:
        raise ValueError("dt must be positive")
    t = np.asarray(np.linspace(0, 120, 6001) if t_eval is None else t_eval, dtype=float)
    if t.ndim != 1 or len(t) == 0 or not np.isfinite(t).all() or np.any(np.diff(t) <= 0) or t[0] < 0:
        raise ValueError("t_eval must be nonempty, finite, nonnegative and strictly increasing")
    for value, name in ((x0, "x0"), (d0, "d0")):
        _finite(value, name)
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must lie in [0,1]")
    initial_a = schedule.coefficients(t[0])[0] if schedule is not None else params.a
    s0 = initial_a if s0 is None else s0
    _nonnegative(s0, "s0")
    current = np.array([x0, s0, d0], dtype=float)
    output = np.empty((len(t), 3))
    output[0] = current
    events = [event for event in schedule.events() if t[0] < event < t[-1]] if schedule else []
    boundaries = [float(t[0]), *events, float(t[-1])] if len(t) > 1 else [float(t[0])]
    work = 0
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        a, r = schedule.coefficients(start + (end - start) / 2) if schedule else (params.a, params.r)
        derivative = lambda y: _rates(model_id, y, params, a, r)
        indices = np.flatnonzero((t > start) & (t <= end))
        if method == "DOP853":
            result = solve_ivp(lambda _, y: derivative(y), (start, end), current,
                               method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True)
            if not result.success:
                raise ArithmeticError(result.message)
            _check_bounds(result.y.T)
            if len(indices):
                output[indices] = result.sol(t[indices]).T
            current = result.y[:, -1].copy()
            work += result.nfev
        else:
            time = start
            destinations = [(float(t[i]), int(i)) for i in indices]
            if not destinations or destinations[-1][0] < end:
                destinations.append((end, None))
            for target, index in destinations:
                while time < target:
                    h = min(dt, target - time)
                    if time + h <= time:
                        raise ArithmeticError("Time step below floating point resolution")
                    current = _rk4_step(current, h, derivative)
                    _check_bounds(current)
                    time = target if h == target - time else time + h
                    work += 4
                if index is not None:
                    output[index] = current
        _check_bounds(current)
    _check_bounds(output)
    return Solution(t.copy(), output, model_id, method, {
        "evidence_class": "synthetic_simulation", "empirical_calibration": False,
        "time_units": "dimensionless_model_time", "state_order": ["x", "s", "d"],
        "active_states": ["x", "s"] if model_id == "M3" else ["x", "d"] if model_id == "M4" else ["x"],
        "s_semantics": "effective_drive_rate; constant placeholder outside M3",
        "events": events, "segments": [[a, b] for a, b in zip(boundaries[:-1], boundaries[1:])],
        "rtol": RTOL if method == "DOP853" else None,
        "atol": ATOL if method == "DOP853" else None,
        "dt_max": dt if method == "RK4" else None, "rhs_evaluations": work,
        "bounds_checked_without_clipping": True,
    })
