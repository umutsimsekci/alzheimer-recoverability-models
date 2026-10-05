# Study design v1: conditional functional recoverability

Date: 2026-10-05. Status: prespecified simulation protocol; **not empirical research or a treatment recommendation**. No ECNP/journal submission. All state, time and intervention scales are dimensionless. Changes to this protocol must be versioned before the affected runs, with a reason; unsuccessful and unresolved cases remain visible.

## 1. Question and estimands

Under which explicit dynamics can a finite repair pulse cause a persistent basin transition after drive reduction, and which measurements distinguish this from monostable relaxation, slow hidden dynamics, structural loss and task-only effects?

Primary estimand: for each fixed synthetic parameter setting, difference in low-impairment target attainment between drive reduction alone and drive reduction plus repair pulse, **under continued drive reduction**. Report raw x and the target outcome separately. The direction of the scalar impairment difference is mathematically imposed; target transition and persistence are the testable quantities.

Secondary estimands: (1) target attainment after **all controls cease**; (2) basin crossing relative to the unstable root at pulse cessation when that root exists; (3) convergence to distinct equilibria from different histories under identical final coefficients; (4) predefined factorial interaction on the x scale, explicitly without a pharmacological-synergy interpretation; (5) cost–attainment frontier in normalized control units; (6) parameter/latent-state distinguishability and held-out prediction under declared observation models. No clinical meaningfulness threshold is asserted.

## 2. Model families

| ID | Dynamics | Purpose |
|---|---|---|
| M0 | `x_dot=a(1-x)-r*x` | Exact linear, monostable comparator. |
| M1 | `x_dot=a(1-x)+b*x^2*(1-x)-r*x` | Legacy bistable candidate; exact root/stability reference. |
| M2 | `x_dot=a(1-x)+b*x*(1-x)-r*x` | Nonlinear monostable comparator; checks whether saturation/nonadditivity is sufficient. |
| M3 | `s_dot=(a_command-s)/tau`, `x_dot=s*(1-x)-r*x` | Slow-input monostable alternative; rate-dependent pseudo-hysteresis. |
| M4 | M1 with `r_effective=r*(1-d)`; `d_dot=epsilon*[x-xc]_+*(1-d)-eta*d` | Structural history; fixed, one-way and reversible variants. |

M0–M3 are primary scalar tests. M4 is a bounded secondary suite. Graph dynamics, PK/PD, patient prediction and optimization of real dosing are deferred. Observation effects are applied separately to these latent dynamics; they are not another unmeasured biological mechanism silently added to improve fit.

## 3. Prespecified design ranges and settings

Core ranges: post-reduction `a in [.03,.09]`, `b in [2,6]`, `r in [.8,1.3]`; pretreatment `a_pre=a+.08`. These are arbitrary legacy-neighborhood ranges chosen for comparison, not estimated parameter intervals. Run the nominal `(a,b,r)=(.06,4,1.05)` and a fixed-seed 256-setting independent uniform design, seed `20261005`. Save the exact parameter table. M0 and M3 ignore b; do not count duplicate b settings as additional biological information. M3 uses `tau in {1,10,50}`, independent of the scalar coefficient sample.

Analytic phase scan: `a in [.001,.16]` with 31 points, `b in [0,6]` with 41 points, `r in [.5,1.8]` with 41 points. Root counts must exclude nonhyperbolic fold cases and handle degeneracies explicitly. This scan illustrates the mathematical region; it is not a parameter likelihood.

Initial-state sensitivity: primary treatment runs start from the largest attracting physical pretreatment equilibrium; label this selection. Two-history tests use `x0=.04` and `.8` at identical fixed coefficients. Repeat treatment tests on the nominal model with `x0 in {.04,.4,.6,.8}`; an exact separatrix start is a numerical stress test, not a stable initial group.

M4 suite: nominal coefficients; `d0 in {0,.15,.30}`, `epsilon in {0,.003,.01}`, `eta in {0,.003,.01}`, `xc=.25`. Record whether the target is initially infeasible or structurally attainable. Separate the `epsilon=eta=0` frozen system from evolving structure.

## 4. Controls and timeline

All groups have the same initial state and treatment onset. Baseline is `t=0..30`; controls begin at `t=30`. Repair pulse adds `.4` to r for `30<=t<45`. Drive reduction sets `a_pre -> a` at 30.

Four arms: no controls; drive reduction alone; repair pulse alone; both. Evaluate each active arm under two policies: maintenance of reduced a through the endpoint (where drive reduction is assigned), or all controls off at 45, including return to `a_pre`. The no-control and repair-only schedules coincide across policies and should be represented once or explicitly marked duplicates.

Primary endpoint `t=120`, persistence interval `[90,120]`. Additional endpoints `240` and `480` are convergence/sensitivity checks; they do not replace the primary result. The nominal full-withdrawal comparison is a planned negative control because `a_pre=.14` is beyond the M1 upper fold. Do not tune the original coefficients to make that case persist.

Cost: `u_a=(a_pre-a_current)/.08` and `u_r=(r_current-r)/.4`, with admissible amplitudes `[0,1]`; abstract cost `C=integral(u_a^2+u_r^2)dt` on the common evaluation horizon. A combination has greater exposure than its components. In the nominal case only, scan repair amplitudes `{0,.2,.4}` and pulse durations `{5,15,30}` under the same policies; report attainable targets versus C. This is an exposure accounting exercise, not toxicity, real affordability or clinical optimal control. Factorial controls and cost reporting do not eliminate the built-in repair ordering.

## 5. Outcomes and targets

For M0–M3, a provisional target is `x<=.15`; for M4 require additionally `d<=.15`. Report `q=(1-x)*(1-d)` as a proposed synthetic readout, not a validated score. Targets must hold at the endpoint and throughout the prescribed persistence interval. Also report sensitivity at x thresholds `.10` and `.20`, and structural thresholds `.10` and `.20` separately. Never pick a threshold after seeing which one favors the pulse.

Continuous outputs: `x(120)`, maximum x during persistence, impairment AUC from 30 to 120, d where applicable, q, maximum drift norm, and distance to any unstable physical equilibrium at cessation. Primary binary target and asymptotic basin classification are distinct. A finite endpoint meeting a target does not prove an attractor transition; a basin transition does not by itself meet the target.

For scalar fixed-coefficient settling, declare convergence only if maximum absolute drift over the last 10 units is `<1e-8` and the nearest residual-certified attracting root is within `1e-6`. Bound the horizon at 480; otherwise report `unresolved_at_horizon`. Root matching alone is insufficient near a fold. For M3 also require its input-state drift. For M4 use the full drift norm, and distinguish equilibrium manifolds from isolated attractors. Negative or null outcomes are ordinary results.

## 6. Hysteresis versus lag experiment

At nominal b,r with frozen structure, ramp commanded a from `.02` to `.14` and back. Use one-leg durations `{60,240,960}`, initial low equilibrium, and relaxation before each run. Handle ramp inputs continuously. Separately compare low/high initial histories at `a=.06` and allow settling to the common convergence limit. Repeat in M0, M2 and M3.

Report loop area and transition locations against ramp speed, but classify bistability from admissible equilibria and history-dependent *settled* endpoints at identical parameters. A loop that contracts with longer relaxation in M3 is a lag control. Slow passage near a genuine fold may also shift transition locations at finite speed, so finite loop width alone is not the discriminator.

## 7. Measurement and identifiability experiment

Measurement designs are explicit synthetic cases:

1. Anchored x only: `y1=1-x+e1`.
2. Redundant tasks: `y2=.85*(1-x)+e2`; no additional latent-state rank is assumed.
3. Anchored functional and structural measurements: `y1=1-x+e1`, `y2=d+e2`.
4. Task-only nuisance: add `.12` to y1 for `30<=t<60`, with no latent change; also test nuisance shared by both tasks.

Use independent Gaussian *observation* errors with SD `{0,.01,.03}`; do not inject unconstrained noise into bounded latent states. Seed `20261006`; 20 noise replicates on the nominal case and 8 fixed design settings selected by index before inspecting results. Sample at 61 equally spaced times from 0 to 120, adding event times if absent. Observation amplitudes are synthetic stress levels.

For anchored scalar x, assess parameter sensitivity rank and singular values, fit-recovery error and profile loss. Demonstrate the steady-state common-rate scaling invariance analytically; do not pretend such data estimate all rates. For M4 explicitly compare known r versus jointly unknown r,d to expose the effective-repair confounding. A known-parameter observability Gramian alone cannot resolve unknown observation maps.

If fitting is implemented within v1, split by complete trajectory/schedule, not individual timepoints: fit baseline/step/pulse records and predict a distinct ramp or withdrawal schedule. Give candidate models the same observation information and document parameter counts/bounds. Compare held-out RMSE in observed units against noise level, residual patterns and parameter profiles. An observationally equivalent or indistinguishable case remains an ambiguity result; no latent restoration claim is accepted solely from fit quality. This fitting substage may be deferred as an explicitly incomplete secondary objective; do not substitute a symbolic rank calculation for performed recovery experiments.

## 8. Numerical accuracy and falsification rules

Implement derivative, equilibria, root residual/stability and exact M0 solution independently from trajectory generation. Adaptive integration with tolerances `rtol=1e-9, atol=1e-11` is the reference; segmented RK4 at `.02,.01,.005` checks event-aligned trajectories. Restart integration exactly at switches; do not round off-grid events silently. Use no clipping. Accepted states satisfy bounds to `1e-9`; root residuals `<1e-10`; away from bifurcations trajectories agree within `1e-6`. Near-fold/separatrix basin labels that disagree under refinement are unresolved, not accepted. Use one-sided analytic handling at the structural threshold.

Concrete falsification/qualification criteria:

- If claimed M1 bistable cases lack three physical simple roots with alternating stable–unstable–stable slopes, reject that bistability classification.
- If nominal maintenance combination does not cross the unstable root and converge near `.070871215` while drive-only converges near `.529128785`, reject reproduction pending numerical diagnosis; do not tune parameters.
- If nominal all-withdrawal retains a low *settled* state at pretreatment `.14`, reject numerical or classification correctness under the declared M1 equations.
- If M0/M2 show two distinct settled physical states at identical final coefficients, reject numerical correctness or the comparator implementation.
- If a claimed hysteresis effect can be reproduced by M3 only as finite-rate lag and vanishes after settling, reject a claim that the observed loop alone identifies bistability.
- If an improvement occurs only in nuisance-contaminated outputs with unchanged x,d, reject latent recovery for that case.
- If sampled synthetic regimes do not show persistent target transitions, report no support in that design; do not enlarge ranges without a new protocol version.
- If unknown-coefficient/observation alternatives remain indistinguishable on held-out schedules, narrow inference to conditional model behavior and state the measurement ambiguity.

These criteria can falsify code/model-specific assertions. They do not by themselves falsify or validate a biological theory of Alzheimer disease.

## 9. Uncertainty, reproducibility and bounded execution

Separate numerical error, design/initial-state dependence, structural model alternatives and observation noise. Uniform-design fractions and bootstrap intervals describe the arbitrary design distribution only; they are not biological uncertainty intervals. Prefer counts with denominators, quantiles and parameter maps; do not attach p-values to deterministic arm differences.

Before a run, save protocol version/hash, exact parameter and control tables, model identifiers, solver settings, package versions, seed, timepoints, code revision where available, and output hashes. Distinguish generated data from reviewed results. Every job records owner, inputs, outputs, dependency, acceptance criteria and elapsed runtime; failures remain explicit. Never overwrite legacy.

Execution order: (A) analytic and core numerical checks; (B) four-arm maintenance/withdrawal and initial-state suite; (C) monostable/ramp/structural alternatives; (D) synthetic observation checks, followed by fitting only if implemented and bounded. Set 10 minutes per initial job and 30 minutes total initial scientific run budget; stop an overrun with a visible incomplete record, profile and split it before retrying. Parallelization may improve throughput but does not increase accepted evidence. No patient experiment, drug selection or clinical intervention is included in v1.
