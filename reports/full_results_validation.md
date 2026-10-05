# Independent validation of full synthetic result artifacts

Date: 2026-10-05. Evidence class: independent artifact/count/hash verification and bounded mathematical checks. No empirical validation or clinical interpretation is supplied. The auditor did not alter drivers, core, tests or result artifacts during this review.

## Current acceptance

**Stage D observation/sensitivity outputs: accepted for the performed scope.** Full observation generation and local sensitivity/algebraic analyses were executed and reviewed. Fit recovery, optimized parameter profiles and candidate-model comparisons on held-out schedules remain unperformed; these are not replaced by local rank or algebraic equivalence curves.

**Full A/B/C outputs: accepted for the performed synthetic scope.** Independent count, hash, endpoint, root, cost and classification checks found no structural or numerical defect in the saved outputs. This is acceptance of reproducible calculations under the declared model and finite schedules, not empirical model validation or completion of the entire prospective protocol.

## Stage D: reviewed artifact directory

`results/observation_full_20261005T083233_0bb911bc`

The manifest says `generated_unreviewed`; this report supplies an external review without mutating that original record. Recomputed SHA256 hashes match **all 12 listed artifacts** and **all 4 source hashes** (protocol, prospective clarifications, observation script and core). No hash mismatch was found.

### Independent coverage and raw-data checks

| Quantity | Independently counted | Expected |
|---|---:|---:|
| Uniform parameter-design rows | 256 | 256 |
| Executed latent cases | 9 | Nominal plus design indices 0..7 |
| Measurement/noise/replicate groups | 2700 | `9*5*3*20` |
| Times per group | 62 | 61 equally spaced times plus off-grid event 45 |
| Raw observation rows | 167400 | `2700*62` |
| Latent trajectory rows | 558 | `9*62` |
| Observation summary rows | 2700 | One per group |
| Local sensitivity records | 36 | `9*4` schedule records |

The saved design agrees exactly with an independent `default_rng(20261005).uniform([.03,2,.8],[.09,6,1.3],(256,3))` reconstruction. Nominal is not counted as an additional random-design setting. Twenty zero-noise repeats are duplicates, and the manifest explicitly marks them as such.

Every raw observation was reconstructed from x, d, saved Gaussian error and the declared nuisance map. Maximum discrepancy was **2.22e-16**. Across all designs/noise levels/repeats, the latent x,d trajectory was reused exactly (maximum discrepancy **0**); noise was not injected into bounded latent states. Recomputed clean-output RMSE, noise-only RMSE, mean first-channel noise and mean nuisance-window offsets match all saved summary rows, maximum discrepancy **5.55e-17**. All groups contain exactly 62 times. The single-channel design correctly omits its second channel in 33480 raw rows.

Nuisance is .12 only for `30<=t<60`. First-channel nuisance and shared raw-amplitude nuisance were distinguished correctly. An observation-only benefit is therefore an explicit mapping example, not a latent intervention effect.

### Independent sensitivity and algebra checks

All eight saved instantaneous observation-Jacobian ranks match independently recomputed matrix ranks. These are known-map instantaneous ranks, not full dynamic observability or global parameter identifiability. The raw equal-amplitude shared-nuisance case has rank two in `(x,d,n)` with d unobserved; the loading-scaled nuisance alternative has rank one. Those cases must not be conflated.

The 36 local trajectory-sensitivity records contain:

- Baseline: 9 rank-one cases.
- Step: 7 rank-three and **2 rank-two** cases.
- Combination: 9 rank-three cases.
- Stacked baseline/step/combination: 9 rank-three cases.

The two rank-two step cases are retained; no universal identifiability claim is justified. Rank uses the declared relative singular-value cutoff 1e-6 and absolute floor 1e-10. The largest reported augmented-state/reference x discrepancy is **5.17e-8**, below the script's 1e-7 criterion.

For the nominal stationary baseline, an independent exact variational solution was used: for scalar slope `f_x` and log-parameter forcing vector `v`, observation sensitivity is `-v*(exp(f_x*t)-1)/f_x`. Its resolved singular value agrees with the saved value to **6.24e-12 relative error**; the exact matrix has rank one. For the nominal step, independently coded event-segmented scalar finite differences at log perturbation 1e-5 agree with all three saved singular values to **4.75e-7 relative error**. The nominal combination finite-difference comparison gives **9.90e-8 relative error**. These checks use fixed known initial x and fixed additive control offsets, matching the declared sensitivity estimand.

Finite differences of an insufficiently resolved near-stationary trajectory can introduce artificial small singular values. The exact stationary variational check, rather than noisy numerical differentiation, supports the baseline rank conclusion. None of these checks performs parameter fitting or establishes practical recovery under observation noise.

Algebraic results are scoped correctly: common rate scaling preserves fixed-input equilibria while changing known-clock transients; fixed-d observations can confound r and d through `r*(1-d)`; direct d measurement or a known additive repair pulse can break that restricted symmetry. The saved r/d equivalence curve is explicitly not an optimized fitted profile. The task-only example retains identical latent counterfactual trajectories and adds .12 to the output; its small numerical state drift is not a biological effect.

## Limits carried forward

These are synthetic design outcomes. Uniform sampling fractions, local ranks, noise simulations, docking tools and LLM agreement are not clinical efficacy evidence. The original study-driver audit's statement that Stage D had not yet run is historical; this later review verifies the performed observation/sensitivity substage while preserving its fitting omissions. The completed A/B/C review follows. The original manifests retain their generation-time statuses; this report supplies review evidence without rewriting those records.

## Full A/B/C: reviewed artifact directory

`results/study_full_v1`

The summary reports 264.621 seconds elapsed and zero recorded failures. Recomputed SHA256 hashes match **all 45 listed output artifacts** and **all four source hashes**. The frozen driver hash is `c11d8991e41b782902cf8a403dc25ffeb01434ef3ec1183b7224bd3d2f32f1aa`. These match the numerical audit's final source. Existing bounded checks independently covered closed-form M0/M2, M3 slow-state responses, M4 regeneration, event segmentation, and separate per-case solves versus the stacked solver; this artifact review does not rerun the entire factorial study.

### Independent counts and phase classification

| Output | Actual | Design interpretation |
|---|---:|---|
| Parameter rows | 257 | 256 uniform settings plus separately labelled nominal |
| Phase grid | 52111 | `31*41*41` |
| Bistable phase points | 5460 | Independent extrema classification agrees at every grid point |
| Numerical check rows | 157 | All pass; largest saved discrepancy 6.35e-8 |
| Primary outcomes | 9252 | `257*6*6`; M0/M1/M2 and M3 tau=1,10,50 |
| Initial-state sensitivity | 144 | Separate sensitivity suite |
| Onset-state sensitivity | 144 | Separate sensitivity suite |
| Factorial contrasts | 3084 | Two policy comparisons per model/tau/setting |
| Structural outcomes | 162 | M4 secondary design |
| Structural extension rows | 288 | Additional horizons, not new independent settings |
| Amplitude/duration/cost outcomes | 18 | Finite pulse family |
| Two-history endpoints | 12 | Two initial states per model/tau |
| Ramp summaries | 18 | Three speeds per model/tau |

The 256 random settings reproduce the seeded uniform design exactly and must not be pooled with the nominal row when reporting design fractions. Each of the 1542 primary `(model,tau,setting)` groups contains six unique schedules: none and repair are each stored once with `identical_policies`, while drive and combination have maintenance and full-withdrawal variants. This avoids counting duplicate no-drive policies as extra evidence.

The independent phase test used scalar derivative extrema, not the driver's root classifier: where `b>3(a+r)`, extrema occur at `(b +/- sqrt(b*(b-3(a+r))))/(3b)` and three physical roots require negative lower-extremum and positive upper-extremum values. There are **zero classification mismatches**, no extrema within 1e-10 of a fold on this grid, and all 1271 b=0 rows are monostable. This grid establishes the model's parameter geometry, not population prevalence.

### Primary endpoints, convergence and imposed ordering

All saved primary endpoint x values are in [0,1] and d remains exactly zero. Independently reconstructed attracting-root distances agree to **2.22e-16**; M1 pre/post bistability flags have **zero mismatches**. Recomputed endpoint vector-field magnitudes are no larger than the recorded maximum over the final ten units, allowing 1e-12 rounding tolerance. Every row labelled converged satisfies the declared saved drift/root-distance rule.

| Model/tau | Converged at 120 | Unresolved at 120 |
|---|---:|---:|
| M0/10 | 1542 | 0 |
| M1/10 | 1481 | 61 |
| M2/10 | 1542 | 0 |
| M3/1 | 1542 | 0 |
| M3/10 | 514 | 1028 |
| M3/50 | 514 | 1028 |

There are 3146 primary extension records. Selecting the latest available record per case yields **8224 converged and 1028 unresolved**, with every remaining unresolved case at horizon 480 in M3/tau=50 under the four schedules that change the drive. The exact slow-state decay is proportional to `exp(-t/50)`; small x drift does not certify convergence while s is still changing. These are finite-horizon unresolved cases, not failures or irreversibility results. Primary target status at 120 and last-horizon convergence answer different questions and must not be substituted for one another.

The minimum paired `x_drive-x_combination` endpoint difference is **-1.33e-15**; no ordering violation exceeds the 1e-6 tolerance. This is expected from the comparison principle when repair is added to the stipulated scalar field. It is not independent evidence of synergy, an effect size in patients, or efficacy of a repair drug.

Among the **256 random M1 settings**, pre/post regime counts are: monostable-to-bistable 84, monostable-to-monostable 170, bistable-to-monostable 2, and bistable-to-bistable 0. The maintenance comparison has **73/256** settings with differing drive-only/combination primary target status, whereas full withdrawal has **0/256**. All other primary model/tau policy comparisons have zero differing target statuses. These fractions describe the specified uniform design and target/window, not uncertainty over biological parameter values.

For nominal M1, the pretreatment and full-withdrawal attracting endpoint is approximately .6027819003. Maintenance drive-only ends at .5291287848; maintenance combination ends at .0708712153 and meets the target. Both full-withdrawal arms end at .6027819003 and fail the target. The post-reduction regime is bistable but the original regime is monostable: return after withdrawal is therefore consistent with the fixed-input geometry. None/repair also return to the same high endpoint.

### Structural state and scope of certification

All 162 structural outcomes and 288 extension rows retain `frozen_d_only` equilibrium scope. Independent endpoint dx/dd calculations agree with the drift bounds; x,d remain in [0,1]. No evolving-structure case is labelled an isolated certified full-system equilibrium. At the latest available horizon across 162 cases:

- 18 are `converged_frozen_structure_manifold`.
- 15 have `full_drift_small_equilibrium_certification_deferred`.
- 129 remain `unresolved_at_horizon`.

Six saved structural rows have borderline .15 threshold flags. Preserve those flags and the drift/certification distinction. Fixed-d roots determine conditional x geometry; they cannot alone certify evolving-d equilibria, establish irreversible loss, or classify a whole structural-control reachable set. Baseline d meeting a target is an initial-state property, not reachability.

### Costs and discriminating alternatives

Every primary cost was independently reconstructed with actual cessation, withdrawal truncation and declared normalization; discrepancy is zero. The same calculation across all 18 pulse-family rows reproduces both total cost and actual duration exactly. Nominal primary costs are 0 (none), 15 (repair), 90 (drive maintenance), 15 (drive withdrawal), 105 (combination maintenance), and 30 (combination withdrawal).

In the secondary family, a planned pulse ending at 60 is truncated at withdrawal 45: actual repair window is 15, not 30. At nonzero amplitudes .2 and .4, all three tested maintenance durations meet the nominal target and all withdrawal variants fail it; amplitude-zero variants fail. The least recorded target-meeting cost in this finite family is 91.25 at amplitude .2 and duration 5, while the nominal .4/15 schedule costs 105. This is not a global minimum over control functions or pulse durations. `actual_repair_duration` denotes the available command window even when amplitude is zero; describe zero-amplitude cases as zero repair exposure.

Constant-input two-history tests distinguish the alternatives: nominal M1 ends at .0708712153 from x0=.04 and .5291287847 from x0=.8, on opposite sides of unstable root .4. M0 and M3 histories converge to .0540540541; M2 histories converge to .7426966706. M3 s is initially equal to the fixed input in this suite, so this particular comparison tests equilibrium history, not slow switching.

Across ramp leg durations 60,240,960, loop areas contract for all models; the M1 areas remain .0492769,.0358320,.0290507, while linear/quadratic/slow-drive loop areas shrink toward zero. Operational x=.4 upward M1 crossing commands are .1335,.115875,.10815625; downward values are .0205,.035625,.0423125. These finite-speed crossing commands are not analytic fold locations. Delayed passage shifts them even in M1, and lag alone generates finite-speed loops in monostable alternatives. M2's initial state already exceeds .4 and is correctly flagged, rather than counted as a new crossing.

## Review disposition and remaining work

**No blocking structural or numerical defect was found in the full saved A/B/C or performed Stage D outputs.** Keep the following limits in any integrated account:

1. Preserve 1028 unresolved primary cases at 480 and 129 unresolved latest-horizon M4 cases; do not relabel them as convergence, failure to recover, or permanent loss.
2. Maintenance target attainment does not establish persistence after all-treatment withdrawal; nominal withdrawal and the full sampled withdrawal comparison are null for a combination target advantage.
3. Extra repair imposes the observed scalar endpoint ordering. The study tests conditional geometry, schedule and observation consequences rather than establishing biological treatment efficacy or synergy.
4. M4 full-system isolated/manifold certification, continuous extrema certification of persistence, global control optimization/reachability, parameter fitting, optimized profiles and held-out candidate-model comparisons remain deferred.
5. Stage D observation generation and local sensitivity/algebraic analyses have completed separately. The A/B/C manifest's broad Stage D deferred entry describes that driver's ownership, not the current integrated completion state. Retain `protocol_complete=false` for the unresolved portions.

All times, rates, pulse amplitudes and costs are dimensionless. Results are uncalibrated synthetic simulations; local sensitivity ranks, synthetic task improvements and sampled fractions do not establish clinical efficacy, safety or realistic biological parameters.
