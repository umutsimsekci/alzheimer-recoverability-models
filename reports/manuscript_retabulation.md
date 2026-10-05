# Saved-output re-tabulation for final manuscript adjudication

Date: 2026-10-05. Evidence: saved synthetic outputs only. No trajectories were rerun, controls changed, thresholds redefined, or primary outcomes replaced. Sources: `results/study_full_v1/primary_outcomes.csv`, `M1_tau10_convergence.csv`, `structural_outcomes.csv`, `structural_convergence.csv`; critique: `logs/critic_manuscript_v1/1/final.txt`.

## Baseline attainment and net gains

All counts below use the 256 random M1 settings, excluding nominal, and the existing sampled x<=.15 window [90,120].

| Schedule | Attainment | Newly attaining relative to no control | Baseline attainers lost |
|---|---:|---:|---:|
| No control | 47/256 | 0 | 0 |
| Maintained drive reduction | 126/256 | 79/256 | 0 |
| Maintained combination | 199/256 | 152/256 | 0 |
| Drive, full withdrawal | 47/256 | 0/256 | 0 |
| Combination, full withdrawal | 47/256 | 0/256 | 0 |

The no-control, drive-withdrawal and combination-withdrawal attaining sets are **identical**, not merely equal-sized. These 47 settings are a subset of both maintenance attaining sets. The primary within-setting combination-versus-drive increment remains 73/256; it is distinct from the gains over no control. Thus the critique's proposed upper bounds of 79 and 152 are exact saved-design counts. Use “new attainment relative to no control within this synthetic design,” rather than biological treatment-attributable benefit.

Suggested addition to abstract/results: “No control already met the target in 47 settings; maintained drive reduction and combination added 79 and 152 attainments, respectively. The combination-versus-drive increment was 73/256.” Do not characterize all 126 or 199 attainments as newly created by intervention.

## Regime membership and finite-window versus settled endpoint

The **73** primary increments comprise **71 monostable-pretreatment to bistable-post-reduction settings**, and **2 monostable-to-monostable settings**. No increment comes from a pretreatment bistable setting. The critique's statement that every increment must involve crossing a post-reduction separatrix is therefore false for the prescribed finite-window estimand.

At t=120, 61 M1 primary rows are unresolved, all random-design rows, representing 28 distinct settings:

| Schedule | Unresolved at 120 |
|---|---:|
| Repair alone | 16 |
| Drive, full withdrawal | 17 |
| Combination, full withdrawal | 17 |
| Maintained drive | 4 |
| Maintained combination | 7 |
| No control | 0 |

All 61 have a later saved record labelled converged by 480. Among these rows, exactly **two** disagree between primary-window attainment and the threshold applied to the latest converged endpoint:

| Setting / arm | Primary [90,120] max x | Primary target | Latest converged endpoint x | Latest time | Endpoint <=.15 |
|---|---:|---|---:|---:|---|
| sample_154 / maintained combination | .1494226921 | Yes | .6947223349 | 480 | No |
| sample_222 / maintained drive | .4732982889 | No | .0582971619 | 240 | Yes |

Both settings have one physical attracting post-reduction root and are the two monostable members of the 73 primary increments. In sample_154, the combination supplies transient window attainment before returning to the high root. In sample_222, repair accelerates attainment of a low root that drive alone reaches later. Neither requires a post-reduction separatrix.

Six of the 73 incremental settings have an unresolved maintenance arm at 120: sample_042 (combination), sample_049 (combination), sample_117 (drive), sample_154 (combination), sample_174 (combination), and sample_222 (drive). The other four retain the same threshold classification at the latest converged endpoint.

For a **separate descriptive endpoint comparison**, use each arm's latest saved converged endpoint: maintained drive has 127/256 below .15 and combination 198/256, an increment of **71/256**. Only sample_154 and sample_222 cease to be incremental; no new incremental setting appears. This endpoint comparison is not the prespecified primary window and must not replace 73/256. It is also not a continuously certified future-persistence result.

Suggested addition: “Of the 73 primary-window increments, 71 occurred in post-reduction bistable settings. Two monostable settings reflected finite-horizon transients; their incremental status disappeared when comparing the latest converged endpoints. The corresponding endpoint increment was 71/256, a different estimand.”

## Frozen and evolving structural counts

“Frozen” means epsilon=eta=0. “Evolving” means at least one is nonzero. The existing maintained-combination joint target is x<=.15 and d<=.15 in the primary window.

| Maintained-combination subset | Target attainment | Borderline cases / attaining | Excluding borderline cases |
|---|---:|---:|---:|
| All structural settings | 16/27 | 1 / 1 | 15/26 |
| Frozen | 2/3 | 1 / 1 | 1/2 |
| Evolving | 14/24 | 0 / 0 | 14/24 |

Borderline here uses the saved x or d .15 borderline flags. The sole borderline maintained-combination row is frozen d0=.15. Six borderline rows exist across the full structural table because the same frozen d0=.15 setting appears under all six schedules. Excluding them is a transparent descriptive sensitivity tabulation, not a revision to target membership. In frozen cases, structural adequacy is an initial-condition filter; in d0=.30, eta=0 cases, impossibility follows directly from d being nondecreasing and the declared joint target.

Across all six schedules there are **162 structural cases: 18 frozen and 144 evolving**. At each case's latest available saved horizon:

- All 18 frozen cases are `converged_frozen_structure_manifold`.
- Of 144 evolving cases, 15 are `full_drift_small_equilibrium_certification_deferred` and 129 are `unresolved_at_horizon`.
- Thus unresolved is **129/162 overall**, or **129/144 evolving**, not “129 of 162 evolving.” Small drift is not certification of an isolated or attracting full-system equilibrium.

Suggested replacement: “Maintained combination met the joint target in 2/3 frozen settings and 14/24 evolving settings. One frozen attainment lay exactly on the d=.15 boundary; excluding flagged cases gives 1/2 frozen and 14/24 evolving. At the last available horizon, 129 of 144 evolving cases remained unresolved; 15 had small full drift with equilibrium certification deferred.”

## Adjudication of analytical shortcuts in the critique

**Reject a universal `exp(-.8*45)` relaxation claim for M1.** A lower bound on r is not a uniform contraction bound. The relevant slope is `f_x=-a-r+2*b*x-3*b*x^2`; it can approach zero near folds or be positive elsewhere. Pretreatment monostability implies the appropriate scalar asymptotic endpoint, but does not by itself establish settling by t=90 or target equivalence throughout a finite window. The observed withdrawal/no-control set equality should be reported as a saved finite-window result explained by final-system geometry, without that unsupported universal rate bound. The two pretreatment bistable settings also require their actual saved outcomes, not inference from monostable counting.

**M3 input-state drift supplies a necessary condition, not a sufficient full-state convergence time.** For maintained reduction, the exact input drift is `(.08/tau)*exp(-(t-30)/tau)`; its 1e-8 crossing times are:

| tau | Maintenance s-drift necessary crossing | Full-withdrawal s-drift necessary crossing |
|---|---:|---:|
| 1 | 45.895 | 60.895 |
| 10 | 165.924 | 178.399 |
| 50 | 629.146 | 576.635 |

For full withdrawal after 45, the drift amplitude is `.08*(1-exp(-15/tau))/tau` and elapsed time is measured from 45. Consequently the critique's approximately 629 time applies only to the maintained-input necessary condition for tau=50, not to every drive-changing schedule. Full convergence additionally requires x drift and attracting-root distance. Because the saved rule takes the maximum over the last ten units, its input-drift component requires an endpoint roughly ten units beyond the listed instantaneous crossing time (for example, beyond 639.146 for maintained tau=50), assuming that full final interval is evaluated. This still does not suffice for full convergence. Do not promise convergence at these algebraic times or extrapolate missing trajectories. Saved records establish that tau=10 cases resolve by 240 and tau=50 active-control cases remain unresolved at 480; the exact input law explains a necessary obstruction.

No biological efficacy, causal mechanism selection, patient probabilities or clinical time scales follow from these re-tabulations.
