"""Stage D synthetic observation audit; no empirical calibration or drug claims."""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
from scipy.integrate import solve_ivp
from recoverability.model import ATOL, RTOL, Parameters, Schedule, equilibria, solve

DESIGN_SEED, OBSERVATION_SEED = 20261005, 20261006
SELECTED_INDICES = tuple(range(8))  # Chosen before inspecting results.
SIGMAS = (0.0, 0.01, 0.03)
MEASUREMENTS = ("anchored_x", "redundant_tasks", "anchored_x_d", "nuisance_y1_only", "nuisance_shared_raw")


def parameter_design():
    """Independent uniform rows; must match the Stage B generator exactly."""
    return np.random.default_rng(DESIGN_SEED).uniform(
        [.03, 2., .8], [.09, 6., 1.3], size=(256, 3))


def observation_times():
    return np.unique(np.r_[np.linspace(0., 120., 61), 30., 45., 60.])


def nuisance(t):
    t = np.asarray(t)
    return .12 * ((t >= 30.) & (t < 60.))


def measurement_values(t, x, d, design):
    if design not in MEASUREMENTS:
        raise ValueError("Unknown measurement design")
    clean = np.column_stack((1. - np.asarray(x), .85 * (1. - np.asarray(x))))
    offset = np.zeros_like(clean)
    if design == "anchored_x_d":
        clean[:, 1] = d
    if design in ("nuisance_y1_only", "nuisance_shared_raw"):
        offset[:, 0] = nuisance(t)
        if design == "nuisance_shared_raw":
            offset[:, 1] = nuisance(t)
    return clean, offset, 1 if design == "anchored_x" else 2


def observation_jacobians():
    matrices = {
        "anchored_x": (["x", "d"], [[-1., 0.]]),
        "redundant_tasks": (["x", "d"], [[-1., 0.], [-.85, 0.]]),
        "anchored_x_d": (["x", "d"], [[-1., 0.], [0., 1.]]),
        "nuisance_y1_only_known_nuisance": (["x", "d"], [[-1., 0.], [-.85, 0.]]),
        "nuisance_shared_raw_known_nuisance": (["x", "d"], [[-1., 0.], [-.85, 0.]]),
        "nuisance_y1_only_unknown_nuisance": (["x", "d", "n"], [[-1., 0., 1.], [-.85, 0., 0.]]),
        "nuisance_shared_raw_unknown_nuisance": (["x", "d", "n"], [[-1., 0., 1.], [-.85, 0., 1.]]),
        "nuisance_shared_loading_scaled_algebraic_alternative": (["x", "d", "n"], [[-1., 0., 1.], [-.85, 0., .85]])}
    records = []
    for name, (columns, matrix) in matrices.items():
        array = np.asarray(matrix)
        records.append({"design": name, "state_columns": columns, "matrix": matrix,
                        "instantaneous_rank": int(np.linalg.matrix_rank(array)),
                        "singular_values": np.linalg.svd(array, compute_uv=False).tolist(),
                        "scope": "instantaneous known-map observation Jacobian, not full dynamic identifiability"})
    return records


def initial_state(params, model="M1", d=0.):
    pre = replace(params, a=params.a + .08)
    stable = [point.x for point in equilibria(model, pre, d=d) if point.stability == "stable"]
    if not stable:
        raise ArithmeticError("No attracting physical pretreatment equilibrium")
    return max(stable)


def make_schedule(params, kind="combination"):
    if kind not in ("baseline", "step", "combination"):
        raise ValueError("Unknown schedule")
    return Schedule(a_pre=params.a + .08, a_post=params.a, r_base=params.r,
                    reduce_drive=kind != "baseline", repair_amplitude=.4 if kind == "combination" else 0.)


def scalar_sensitivities(params, schedule, x0, times):
    """M1 sensitivities of y=1-x wrt log(a,b,r), fixed known x0.

    The additive .08 pre-drive offset and .4 repair pulse are fixed known commands.
    Their amplitudes are not scaled by the perturbed baseline parameters.
    """
    state = np.r_[x0, 0., 0., 0.]
    output = np.empty((len(times), 4))
    output[0] = state
    events = [v for v in schedule.events() if times[0] < v < times[-1]]
    boundaries = [float(times[0]), *events, float(times[-1])]
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        a_current, r_current = schedule.coefficients((left + right) / 2.)
        def derivative(_, z):
            x, sensitivities = z[0], z[1:]
            drift = a_current * (1. - x) + params.b * x*x * (1. - x) - r_current*x
            fx = -a_current + params.b * (2.*x - 3.*x*x) - r_current
            forcing = np.array([params.a * (1. - x), params.b*x*x*(1. - x), -params.r*x])
            return np.r_[drift, fx*sensitivities + forcing]
        result = solve_ivp(derivative, (left, right), state, method="DOP853",
                           rtol=RTOL, atol=ATOL, dense_output=True)
        if not result.success:
            raise ArithmeticError(result.message)
        indices = np.flatnonzero((times > left) & (times <= right))
        output[indices] = result.sol(times[indices]).T
        state = result.y[:, -1]
    reference = solve("M1", params, schedule, x0=x0, t_eval=times)
    error = float(np.max(np.abs(output[:, 0] - reference.x)))
    if error > 1e-7:
        raise ArithmeticError(f"Augmented/reference latent disagreement: {error}")
    return -output[:, 1:], error


def rank_record(matrix):
    values = np.linalg.svd(matrix, compute_uv=False)
    relative_cutoff = 1e-6
    threshold = max(1e-10, relative_cutoff * values[0])
    return {"rank": int(np.count_nonzero(values > threshold)), "singular_values": values.tolist(),
            "relative_cutoff": relative_cutoff, "absolute_floor": 1e-10,
            "condition_number": float(values[0] / values[-1]) if values[-1] > threshold else None,
            "columns": ["log_a", "log_b", "log_r"],
            "scope": "local known-map M1 trajectory sensitivity; no global identifiability or performed fit"}


def algebraic_examples(times):
    p = Parameters(.06, 4., 1.05)
    d_a, d_b = .15, .25
    effective = p.r * (1. - d_a)
    alternate = replace(p, r=effective / (1. - d_b))
    first = solve("M4", p, x0=.4, d0=d_a, t_eval=times)
    equivalent = solve("M4", alternate, x0=.4, d0=d_b, t_eval=times)
    known_r_alt = solve("M4", p, x0=.4, d0=d_b, t_eval=times)
    profile = []
    for d in (0., .1, .15, .2, .3, .4):
        r = effective / (1. - d)
        other = solve("M4", replace(p, r=r), x0=.4, d0=d, t_eval=times)
        loss = float(np.mean((other.x - first.x)**2))
        profile.append({"d": d, "r": r, "effective_r": r*(1.-d), "x_mean_squared_difference": loss})
    # A fixed known additive pulse acts on r before multiplying by (1-d): it breaks this symmetry.
    def pulse(par):
        return Schedule(par.a, par.a, par.r, reduce_drive=False, repair_amplitude=.4)
    pulse_a = solve("M4", p, pulse(p), x0=.4, d0=d_a, t_eval=times)
    pulse_b = solve("M4", alternate, pulse(alternate), x0=.4, d0=d_b, t_eval=times)
    roots_a = [v.x for v in equilibria("M1", p)]
    scaled = Parameters(2.*p.a, 2.*p.b, 2.*p.r)
    roots_b = [v.x for v in equilibria("M1", scaled)]
    dynamic_a = solve("M1", p, x0=.2, t_eval=times)
    dynamic_b = solve("M1", scaled, x0=.2, t_eval=times)
    steady = {"scale": 2., "identity": "F(x;c*a,c*b,c*r)=c*F(x;a,b,r)",
              "equilibrium_roots_original": roots_a, "equilibrium_roots_scaled": roots_b,
              "max_root_difference": float(np.max(np.abs(np.array(roots_a)-roots_b))),
              "max_clock_time_trajectory_difference": float(np.max(np.abs(dynamic_a.x-dynamic_b.x))),
              "interpretation": "Fixed-input steady-state data cannot identify the common rate scale. Known-clock transient data can distinguish it; fixed additive intervention amplitudes also break simple scaling."}
    # Pure observation change: no control, constant frozen latent equilibrium throughout.
    stationary_x = initial_state(p, "M4", d_a)
    stationary = solve("M4", replace(p, a=p.a+.08), x0=stationary_x, d0=d_a, t_eval=times)
    clean, offset, _ = measurement_values(times, stationary.x, stationary.d, "nuisance_shared_raw")
    task_only = {"x0": stationary_x, "latent_max_drift_from_initial": float(np.max(np.abs(stationary.x-stationary_x))),
                 "d_max_drift_from_initial": float(np.max(np.abs(stationary.d-d_a))),
                 "observed_nuisance_amplitude": .12, "window": "30<=t<60",
                 "latent_counterfactual_difference": 0.,
                 "interpretation": "Same latent trajectory, different observation map; apparent task improvement is not latent recovery."}
    confounding = {"r_a": p.r, "d_a": d_a, "r_b": alternate.r, "d_b": d_b,
                   "effective_r": effective, "equivalence_formula": "r_b=r_a*(1-d_a)/(1-d_b)",
                   "unknown_r_d_jacobian": [[1.-d_a, -p.r]], "rank": 1,
                   "known_r_d_derivative": -p.r,
                   "max_anchored_x_difference_equivalent_pair": float(np.max(np.abs(first.x-equivalent.x))),
                   "max_anchored_d_difference_equivalent_pair": abs(d_a-d_b),
                   "max_x_difference_known_r_changed_d": float(np.max(np.abs(first.x-known_r_alt.x))),
                   "max_x_difference_known_additive_pulse_pair": float(np.max(np.abs(pulse_a.x-pulse_b.x))),
                   "algebraic_loss_curve_not_fitted_profile": profile,
                   "scope": "Frozen structure, constant coefficients or identical effective repair. Other coefficients and x0 known. No empirical estimation."}
    ambiguity = {"observation_map":"[1-x+n,.85*(1-x+n)]", "pair_a":{"x":.5,"n":0.},
                 "pair_b":{"x":.6,"n":.1}, "identical_observation_pair":[.5,.425],
                 "scope":"Algebraic loading-scaled nuisance alternative; not the primary identical raw-amplitude shared-nuisance simulation."}
    return {"steady_state_scaling": steady, "frozen_m4_confounding": confounding,
            "task_only_false_improvement": task_only, "loading_scaled_nuisance_ambiguity":ambiguity}, (stationary, clean, offset), (first, equivalent, known_r_alt)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", "utf-8")


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_figures(directory, examples):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    algebra, stationary_packet, pairs = examples
    stationary, clean, offset = stationary_packet
    figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(stationary.t, stationary.x, label="latent x, unchanged")
    axes[0].plot(stationary.t, stationary.d, label="frozen d, unchanged")
    axes[0].set_ylabel("Synthetic latent units"); axes[0].legend()
    axes[1].plot(stationary.t, clean[:, 0], label="clean y1")
    axes[1].step(stationary.t, clean[:, 0]+offset[:, 0], where="post", label="y1 with task nuisance")
    axes[1].axvspan(30, 60, color="orange", alpha=.12)
    axes[1].set_xlabel("Dimensionless model time"); axes[1].set_ylabel("Synthetic observed units"); axes[1].legend()
    figure.suptitle("Observation-only task improvement; no latent change")
    figure.tight_layout(); figure.savefig(directory / "task_only_false_improvement.png", dpi=180); plt.close(figure)
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    for trajectory, label in zip(pairs, ("r=1.05, d=.15", "unknown r adjusted, d=.25", "known r=1.05, d=.25")):
        axes[0].plot(trajectory.t, 1.-trajectory.x, label=label)
    axes[0].set_xlabel("Dimensionless model time"); axes[0].set_ylabel("Anchored 1-x"); axes[0].legend(fontsize=8)
    profile = algebra["frozen_m4_confounding"]["algebraic_loss_curve_not_fitted_profile"]
    axes[1].plot([v["d"] for v in profile], [v["r"] for v in profile], "o-")
    axes[1].set_xlabel("Frozen unknown d"); axes[1].set_ylabel("Equivalent r=k/(1-d)")
    axes[1].set_title("Exact effective-repair equivalence curve")
    figure.suptitle("Synthetic frozen-M4 ambiguity; no fitted recovery experiment")
    figure.tight_layout(); figure.savefig(directory / "frozen_m4_confounding.png", dpi=180); plt.close(figure)


def run(mode="quick", output_root=None, timeout_seconds=600):
    if mode not in ("quick", "full") or not 1 <= timeout_seconds <= 600:
        raise ValueError("mode quick/full; timeout must be 1..600 seconds")
    started = time.monotonic()
    output_root = Path(output_root) if output_root else ROOT / "results"
    directory = output_root / ("observation_" + mode + "_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "_" + uuid.uuid4().hex[:8])
    directory.mkdir(parents=True, exist_ok=False)
    times, design = observation_times(), parameter_design()
    cases = [("nominal", -1, Parameters(.06, 4., 1.05))] + [
        (f"design_{i:03d}", i, Parameters(*design[i])) for i in SELECTED_INDICES]
    executed_cases = cases if mode == "full" else cases[:2]
    repetitions = 20 if mode == "full" else 2
    expected_full = len(cases)*len(MEASUREMENTS)*len(SIGMAS)*20
    expected_mode = len(executed_cases)*len(MEASUREMENTS)*len(SIGMAS)*repetitions
    manifest = {"evidence_class": "synthetic_simulation", "empirical_calibration": False,
                "mode": mode, "status": "running", "design_seed": DESIGN_SEED,
                "observation_seed": OBSERVATION_SEED, "prespecified_design_indices": list(SELECTED_INDICES),
                "times": times.tolist(), "time_units": "dimensionless_model_time",
                "noise": {"distribution": "independent Gaussian observation errors", "sd": list(SIGMAS),
                          "full_replicates": 20, "executed_replicates": repetitions,
                          "rng_key": "SeedSequence([20261006, design_index+1, measurement_index, sd_index, replicate])",
                          "noiseless_replicates": "duplicates, not additional independent evidence"},
                "latent_observation_suite": {"model": "M4", "d0": .15, "epsilon": 0., "eta": 0.,
                                             "schedule": "maintenance combination", "solver": "DOP853", "rtol": RTOL, "atol": ATOL},
                "source_hashes": {"protocol": hash_file(ROOT/"protocol/study_design_v1.md"),
                                  "prospective_clarifications_v1_1": hash_file(ROOT/"protocol/clarifications_v1_1.md"),
                                  "script": hash_file(Path(__file__)), "core": hash_file(ROOT/"src/recoverability/model.py")},
                "versions": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")},
                "fitting_substage": {"status": "deferred", "reason": "Full fit recovery, optimized profiles and schedule-held-out candidate comparisons are not implemented within this bounded Stage D audit. Sensitivities and algebraic equivalence curves are not substitutes."}}
    write_json(directory / "manifest.json", manifest)
    with (directory / "parameter_design_256.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream); writer.writerow(["index", "a", "b", "r"])
        writer.writerows((i, *row) for i, row in enumerate(design))
    write_json(directory / "selected_cases.json", [{"case": c, "index": i, "a": p.a, "b": p.b, "r": p.r} for c,i,p in cases])
    write_json(directory / "observation_jacobians.json", observation_jacobians())
    completed, actual_observation_rows, failures, sensitivity_records, summary_rows = 0, 0, [], [], []
    fields = ["case", "design", "sigma", "replicate", "time", "x", "d", "y1", "y2", "e1", "e2", "nuisance1", "nuisance2"]
    controls = []
    try:
        with (directory / "observations.csv").open("w", newline="", encoding="utf-8") as stream, (directory / "latent_trajectories.csv").open("w", newline="", encoding="utf-8") as latent_file:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
            latent_writer = csv.writer(latent_file); latent_writer.writerow(["case", "time", "x", "d"])
            for case, index, params in executed_cases:
                if time.monotonic() - started > timeout_seconds - 15:
                    failures.append({"case": case, "reason": "bounded runtime reached before case"}); break
                try:
                    schedule = make_schedule(params)
                    x0 = initial_state(params, "M4", .15)
                    latent = solve("M4", params, schedule, x0=x0, d0=.15, t_eval=times)
                    controls.append({"case": case, "a_pre": params.a+.08, "a_post": params.a, "r_base": params.r,
                                     "onset": 30., "pulse_end": 45., "repair_amplitude": .4, "withdrawal": None,
                                     "x0": x0, "initial_selection": "largest attracting pretreatment equilibrium, frozen d=.15"})
                    latent_writer.writerows((case,t,x,d) for t,x,d in zip(times,latent.x,latent.d))
                    for design_index, measurement in enumerate(MEASUREMENTS):
                        clean, offset, channels = measurement_values(times, latent.x, latent.d, measurement)
                        for sigma_index, sigma in enumerate(SIGMAS):
                            for replicate in range(repetitions):
                                rng = np.random.default_rng(np.random.SeedSequence([OBSERVATION_SEED,index+1,design_index,sigma_index,replicate]))
                                errors = rng.normal(0., sigma, size=clean.shape)
                                observed = clean + offset + errors
                                for j,t in enumerate(times):
                                    writer.writerow(dict(zip(fields, [case,measurement,sigma,replicate,t,latent.x[j],latent.d[j],observed[j,0],
                                        observed[j,1] if channels==2 else "", errors[j,0],errors[j,1] if channels==2 else "",offset[j,0],offset[j,1] if channels==2 else ""])))
                                    actual_observation_rows += 1
                                mask = (times>=30.) & (times<60.)
                                summary_rows.append({"case":case,"design":measurement,"sigma":sigma,"replicate":replicate,
                                    "observation_rmse_vs_clean": float(np.sqrt(np.mean((observed[:,:channels]-clean[:,:channels])**2))),
                                    "noise_rmse_after_nuisance_removed": float(np.sqrt(np.mean(errors[:,:channels]**2))),
                                    "mean_y1_noise":float(np.mean(errors[:,0])), "window_mean_y1_offset":float(np.mean(offset[mask,0])),
                                    "latent_nuisance_counterfactual_max_difference":0.})
                                completed += 1
                    fixed_x0 = initial_state(params)
                    matrices = []
                    for kind in ("baseline", "step", "combination"):
                        matrix, error = scalar_sensitivities(params, make_schedule(params,kind), fixed_x0,times)
                        matrices.append(matrix)
                        sensitivity_records.append({"case":case,"schedule":kind,"fixed_known_x0":fixed_x0,
                                                    "reference_max_x_error":error, **rank_record(matrix)})
                    sensitivity_records.append({"case":case,"schedule":"stacked_baseline_step_combination",
                                                "fixed_known_x0":fixed_x0,**rank_record(np.vstack(matrices))})
                except Exception as error:
                    failures.append({"case":case,"reason":f"{type(error).__name__}: {error}"})
        examples = algebraic_examples(times)
        write_json(directory / "algebraic_examples.json", examples[0])
        write_figures(directory, examples)
    except Exception as error:
        failures.append({"stage":"output_or_algebra","reason":f"{type(error).__name__}: {error}"})
    write_json(directory / "controls.json", controls)
    write_json(directory / "parameter_sensitivity.json", sensitivity_records)
    with (directory / "observation_summary.csv").open("w",newline="",encoding="utf-8") as stream:
        columns = list(summary_rows[0]) if summary_rows else ["case"]
        writer = csv.DictWriter(stream,fieldnames=columns); writer.writeheader(); writer.writerows(summary_rows)
    coverage = {"expected_full_measurement_replicates":expected_full,"expected_mode_measurement_replicates":expected_mode,
                "actual_measurement_replicates":completed,"expected_mode_observation_rows":expected_mode*len(times),
                "actual_observation_rows":actual_observation_rows,"full_case_count":len(cases),"executed_case_count":len(controls),
                "expected_mode_sensitivity_records":4*len(executed_cases), "actual_sensitivity_records":len(sensitivity_records),
                "full_mode_complete":mode=="full" and completed==expected_full and not failures,
                "mode_complete":completed==expected_mode and not failures,
                "fitting_substage_complete":False,"failures":failures}
    write_json(directory / "expected_vs_actual.json",coverage)
    manifest.update({"status":"generated_unreviewed" if coverage["mode_complete"] else "incomplete_unreviewed",
                     "elapsed_seconds":round(time.monotonic()-started,3),"coverage":coverage,
                     "artifact_hashes":{p.name:hash_file(p) for p in sorted(directory.iterdir()) if p.is_file() and p.name!="manifest.json"}})
    write_json(directory / "manifest.json",manifest)
    return directory, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=("quick","full"),default="quick")
    parser.add_argument("--output-root",type=Path)
    parser.add_argument("--timeout-seconds",type=int,default=600)
    args = parser.parse_args()
    directory, manifest = run(args.mode,args.output_root,args.timeout_seconds)
    print(json.dumps({"result_directory":str(directory),"status":manifest["status"],
                      "elapsed_seconds":manifest["elapsed_seconds"],"coverage":manifest["coverage"]}))
    return 0 if manifest["coverage"]["mode_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
