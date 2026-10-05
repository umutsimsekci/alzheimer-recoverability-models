# Independent experiment-driver validation

Date: 2026-10-05. Scope: `scripts/run_study_v1.py`, protocol fidelity and bounded numerical checks. Evidence class: mathematical/technical verification of synthetic simulations, not empirical validation. The auditor authored the small core earlier; therefore the driver checks include independently coded closed-form references and error-injection tests as well as core comparisons. The driver and core were not modified by this auditor.

## Acceptance status

**Accepted for the bounded mathematical/numerical checks below.** After the correction described below, actual execution of `.venv/Scripts/python.exe -m pytest tests/test_study_driver.py -q` returned **19 passed in 7.76 seconds**, exit code 0. This accepts the driver checks, not a complete protocol result or empirical theory. Stage D remains incomplete. The prior failed check and correction are retained as evidence.

## Defect and corrective requirement

A stacked M1 solve containing the nominal setting and settings immediately on either side of the upper fold (`a=.1027147093197253 +/- 1e-5`) differs from independent per-case solves by **1.44879866e-6**, exceeding the protocol's `1e-6` trajectory comparison criterion. The nominal no-control trajectory exhibits about 1e-6 error around its known equilibrium near t=123 in the mixed batch. Shared RMS error control over many independent state components can conceal a case's error. This is a numerical stress test, not an expanded scientific parameter range.

The driver author corrected this by bounding DOP853 to `max_step=.5` and recording the bound in the manifest. Model coefficients and the acceptance threshold remained unchanged. The subsequent full 19-test validation passed, including the formerly failing case. The standard full 1542-case M1 stack also passed the prespecified selected-case comparisons. This defect/correction concerns solver robustness, not biological evidence.

## Verified numerical comparisons

The checks below were executed; they are not proposed future tests:

| Check | Reference independent of the driver | Final result |
|---|---|---|
| M0 off-grid pulse end/withdrawal | Exact eventwise exponential solution, separately coded | Passed at absolute tolerance 1e-7. |
| M2 off-grid pulse end/withdrawal | Closed-form Riccati solution, separately coded from its two quadratic roots | Passed at absolute tolerance 1e-7. |
| M3 slow commanded drive | Exact exponential step response and exact piecewise linear-ramp response | Passed at absolute tolerance 1e-9. |
| M4 structural regeneration | Exact `d(t)=d0*exp(-eta*t)` for zero impairment/damage generation | Passed at absolute tolerance 1e-9. |
| Full M1 design stack | 1542 cases (257 settings x 6 distinct schedules); separately solved nominal and sample indices 0,1,7,127,255, all 6 arms | Passed max absolute 1e-6 over the complete 0..480 output trajectory, including endpoints 120 and 480. |
| M0/M2/M3 and near-fold M1 batches | Individual core solves with the same event schedule | Passed at max absolute 1e-6 after the bounded-step correction. |
| M4 frozen and evolving mixtures | Individual core solves for fixed structure, one-way loss and regeneration variants | Passed max absolute 1e-6 over 0..480. |
| Internal accepted steps | Injected hidden out-of-bound x/d and negative s while sampled output remains valid | Driver rejected all three injected cases. |

Per-case core comparisons reduce the batch-RMS concern but are not completely independent mathematical validation because the core was authored by this auditor. Analytic M0/M2/M3/M4 special cases supply additional independent references. Neither solver agreement nor an exact solution validates model biology.

## Protocol fidelity and interpretation

- The fixed seed and 256-setting design retain the nominal case separately. The six distinct schedules explicitly merge no-control and repair-only maintenance/withdrawal duplicates under `identical_policies`; factorial contrasts reuse those cases with correct interpretation labels.
- Withdrawal restores drive and ends repair. A planned 30-unit repair pulse with withdrawal at 45 is truncated at 45; it does not prolong all-treatment exposure. Nominal abstract costs over 0..120 are `[0,15,90,15,105,30]` for the six schedules. These values passed a direct test.
- Nominal combination under maintenance ends near x=.070871215; all-treatment withdrawal returns near x=.602781900. The former meets the sampled x target and the latter does not. These are synthetic states, not patient recovery.
- Convergence requires the largest absolute active drift during the sampled last ten time units below 1e-8 and distance to an attracting root below 1e-6. A manufactured trajectory that reaches a root only at the final point is correctly unresolved; an unsettled M3 support state is also unresolved.
- M4 equilibrium flags explicitly use frozen d. Evolving M4 uses the label `full_drift_small_equilibrium_certification_deferred` rather than asserting an isolated asymptotically stable full equilibrium. Frozen M4 labels its neutral structural direction/manifold. Full evolving-structure equilibrium certification remains deferred; unresolved structural cases now receive 240/480 extensions in full mode, without converting small drift into certified full equilibrium stability.
- Two fixed-input histories settle to distinct nominal M1 roots; M0/M2/M3 histories settle together. Ramp output describes rate-dependent loops and threshold crossings, not observed biological folds. M3 ramp support agrees with its analytic forcing solution.
- Target persistence is checked at sampled .25 intervals plus events. Continuous maxima are not certified, and the manifest marks that limitation. Report sampled target attainment explicitly; do not convert it into a mathematical guarantee throughout the entire interval.
- The driver records failed/incomplete jobs, exception text and tracebacks, output hashes and unreviewed status. An injected Stage A failure returns a failure exit code and remains visible; attempts to overwrite its output directory raise FileExistsError.

## Required reporting limits

`--quick` is a reduced technical smoke, not full v1 science. `--full` still carries `protocol_complete=False`: Stage D observation/noise/identifiability/fitting and several secondary certifications are deferred. Complete manuscript claims must preserve those omissions. Root counts, costs, two-history results and terminal target flags are conditional on the stated dynamics and arbitrary parameter design. No clinical efficacy, mechanistic proof or drug-pair synergy is supported. The later metadata-only root-feasibility/margin additions do not expand this acceptance into completed simulation runs; run artifacts still require review.
