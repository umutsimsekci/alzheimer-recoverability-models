# Reading and reproducing the 19 October inference revision

This v2.0.0 reader package adds a synthetic inference benchmark. The v1.0.1 tag, original core equations, original protocols and saved dynamics outputs are unchanged. The new work is synthetic and uncalibrated; it does not identify a biological mechanism in patients or establish treatment efficacy. Code availability does not imply journal acceptance.

## What the new experiment tests

The comparison uses M0–M3 as both data-generating and fitted candidate families. Training design A contains two independent noisy step traces. Design B contains the same first step trace and a step-plus-repair trace; second-trace noise is paired across designs. Each design has 132 observations. Two whole schedules are held out: an earlier, smaller pulse and complete withdrawal. Their observations never select a model, start or parameter bound. This is distinct from withholding random timepoints within a fitted schedule.

Three parameter settings, two known noise levels and four replicates produce 192 datasets and 768 candidate fits. The state starts at a known x(0)=0.6; the observation is the known anchored map y=1−x with known Gaussian noise SD. Input timing and additive amplitudes are also known. These favourable assumptions limit the inference. Family selection is only among the four fitted alternatives, and parameters are dimensionless hypotheses.

The additional frozen-M4 experiment estimates r and fixed d with a=0.06, b=4, initial state and observation map held known. It contrasts two step traces, step plus known repair, and step plus direct d measurements at equal scalar observation counts. It does not estimate evolving damage or unknown input/measurement mechanisms.

## Files and interpretation

All new numerical files are under `revision_20261019/analysis/`.

| Files | Contents and reading rule |
| --- | --- |
| `M0/` … `M3/observations.csv` | Saved training and withheld noisy observations and their clean generators; split labels distinguish their roles. |
| `M0/` … `M3/fits.json`, `fits.csv`, `fit_checkpoints.jsonl` | Every candidate and all deterministic starts, convergence/failure, parameter values, training criterion, boundary diagnostics and prediction-error summaries. JSON preserves nested start records. |
| `M0/` … `M3/heldout_predictions.csv` | Complete withheld trajectories for successful fitted candidates. These predictions were not used for selection. |
| `frozen_M4/` | Observations, joint fit records, checkpoints, optimized grid profiles and the generation manifest. |
| `selection_summary.csv` | All-assigned primary selection counts and selected-model prediction errors. |
| `candidate_failures.csv` | All nine failed candidates, including five generating-family candidates, with optimizer messages. |
| `parameter_recovery_summary.csv` | Relative parameter errors conditional on successful true-family fitting; `assigned_n`, `failed_n` and `n` give the denominator. |
| `paired_sensitivity_summary.csv` | All-assigned paired comparison and explicitly post hoc convergence-conditioned subsets. |
| `post_hoc_M1_paired_prediction.csv` | Post hoc comparison fixing the fitted family to M1 and requiring successful M1 fits in both designs. It cannot replace the primary comparison. |
| `frozen_profile_summary.csv`, `profile_resolution_flags.csv` | Grid-set extents and boundary/resolution qualifications. |
| `numerical_validation.csv`, `numerical_validation.json` | Forward-solver and sensitivity checks, distinct from optimizer success and biological validity. |
| `summary.json` | Combined saved summaries, qualifications, failed-fit and postprocessing records. |
| `figures/` | Three data-derived figures as PDF/PNG/SVG and CSV source values. These PDFs are original figures, not copied publications. |

There were 37 failed starts among 2304 starts. Nine of 768 candidates had no successful start. The failed candidates remain present; they do not receive finite model-selection criteria. Five failures involve the generating family, all under design A. As a result, a design difference can reflect both statistical information and practical optimizer reliability. The post hoc subsets make this limitation visible without changing the prespecified primary dataset or replacing failed fits. Parameter-error summaries based only on successful true-family fits require their denominators.

The held-out comparison is to both known clean generator trajectories and independently noisy observations. A small error tests prediction within the declared synthetic system. It does not establish empirical generalization. Numerical checks compare forward solutions and sensitivities, including independent reference integration; they do not prove global optimization or biological validity.

## Frozen-M4 ambiguity and grid limits

Step-only data depend on r and d through the effective repair r(1−d). The true product is 0.8925. Along the equivalent family r(d)=0.8925/(1−d), the imposed r≤4 bound limits d to 0.776875. Thus the step-only grid endpoint near 0.775 is primarily a parameter-bound consequence, not information identifying d.

The known additive repair pulse can break this particular symmetry, as can a calibrated direct-d channel. The profile computes optimized loss at 33 points with spacing 0.025 and uses a descriptive change-in-standardized-SSE threshold of 3.841. Threshold crossings are not refined and repeated-sampling coverage is not established. A single accepted grid point means the width is unresolved on this grid; it does not imply zero uncertainty. The adjacent enclosing bracket may be 0.05 wide. No conclusion here extends automatically to unknown a or b, evolving d, unknown perturbation amplitudes, or unknown observation mappings.

## Exact executed source versus the reader runner

The immutable file `revision_20261019/analysis/inference_source_executed.py` is the byte-exact source used for the saved scientific fits:

```text
SHA256 6b9f1e54cc1b472957a8255ede914be607596f9987e15414eb502e74f3cf1a84
```

Do not execute that snapshot in its storage directory: its root-relative imports assume the `scripts/` location. The usable `scripts/run_inference_revision_20261019.py` differs only by removing the obsolete output field `diagnostic_noise_free_truth_rmse=0`. That historical zero in the saved records is an identity of the known generating map, not evidence of noiseless parameter recovery. It is excluded from all result summaries and figures. Fitting, selection and numerical methods are otherwise unchanged. A small unit test of noiseless M0 fitting is a technical check, not a full noiseless inference experiment.

The original manifests retain the executed source hash and generation-time `complete_unreviewed` status. The aggregate summary retains `generated_pending_independent_review`. These are historical program states, not claims that a journal reviewed or accepted the work. `docs/inference_source_provenance.json` lists the exact public source copies and hashes, including the separate postprocessor. No private manuscript, journal assessment, account data or third-party paper is included in this addition.

## Reproduction

Use the Python environment and exact package versions in the root README and `requirements.lock.txt`. Original inference execution used Python 3.13.6, NumPy 2.5.3 and SciPy 1.18.1; figure generation also uses Matplotlib. No GPU, API key or provider account is needed. Commands below run from the repository root. On Unix substitute `.venv/bin/python` for `.venv/Scripts/python.exe`.

First run the independent technical tests; they take neither a clinical dataset nor a full grid run:

```text
.venv/Scripts/python.exe -m pytest -q -p no:cacheprovider tests/test_inference_revision_20261019.py
```

For fresh numerical output choose a new folder, preserving the supplied data. The script default points at the saved analysis directory, so the explicit `--output` is intentional. Each fit/profile job has its own time bound; slower machines can reach that bound, which must remain a recorded failure rather than being hidden. The smoke is a runtime/forward check, not the full design.

```text
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode validate --output results/reproduced_inference
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode smoke --generator M1 --output results/reproduced_inference --timeout 600
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode run --generator M0 --output results/reproduced_inference --timeout 600
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode run --generator M1 --output results/reproduced_inference --timeout 600
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode run --generator M2 --output results/reproduced_inference --timeout 600
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode run --generator M3 --output results/reproduced_inference --timeout 600
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode profiles --output results/reproduced_inference --timeout 600
.venv/Scripts/python.exe scripts/run_inference_revision_20261019.py --mode summarize --output results/reproduced_inference
```

Run postprocessing only after all four full generator jobs and frozen profiles complete. Copy the postprocessor into the new output folder because it deliberately reads adjacent results. Save the runner actually used in that folder under the expected snapshot name, so its provenance describes the fresh run rather than falsely reusing the original hash:

```text
.venv/Scripts/python.exe -c "from pathlib import Path; import shutil; p=Path('results/reproduced_inference'); shutil.copy2('scripts/run_inference_revision_20261019.py',p/'inference_source_executed.py'); shutil.copy2('revision_20261019/analysis/summarize_and_plot.py',p/'summarize_and_plot.py')"
.venv/Scripts/python.exe results/reproduced_inference/summarize_and_plot.py
```

The postprocessor performs no new fits or model selection; it labels post hoc diagnostics and draws figures from saved results. For exact reproduction of the historical output schema, work in a disposable copy and replace the runner there with the immutable executed snapshot before starting. Never overwrite the original analysis files or retag v1.0.1. Timestamps, elapsed time and rendered figure bytes may differ across environments even when numerical values agree. Cross-platform execution and a fresh dependency installation have not been verified by this packaging step.

The earlier 1028 unresolved M3 cases and 129 unresolved evolving-M4 cases remain unresolved at their saved horizons. The inference addition supplies no new asymptotic certification, empirical calibration, clinical prediction or treatment recommendation.
