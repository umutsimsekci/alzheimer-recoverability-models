# Reviewed synthetic observation audit

Full run: results/observation_full_20261005T083233_0bb911bc. Date: 5 October 2026. Coordinator review: Codex. Evidence class: synthetic computation and algebra, not empirical estimation.

Eight focused unit tests passed before execution. Quick mode then completed 60/60 packets and 8/8 local sensitivity records. Full mode completed all 2700/2700 measurement packets (167400 observation rows), 9/9 selected cases, and 36/36 sensitivity records, with zero recorded failures. The 20 repeated zero-noise packets per case/design are exact duplicates and not additional independent evidence. Full generation took approximately 4.225 seconds.

The fixed-input r=1.05,d=.15 and r=1.19,d=.25 pair had identical anchored x trajectories numerically, while anchored d differed by .1. With r held known and d changed, maximum x difference was about .09578. Applying the same known additive repair pulse to both r values broke this particular symmetry (maximum x difference about .08106). No global identifiability conclusion follows from those examples.

Doubling all autonomous scalar rates preserved the equilibrium roots exactly while changing a known-clock transient. In the nominal fixed-known-initial-state sensitivity calculation, baseline local rank was 1; step and combination local ranks were 3 under the declared cutoff. Their condition numbers were approximately 4865 and 195 respectively; rank alone does not ensure practical parameter recovery under noise. The combined schedule stack had condition number approximately 113. These are local sensitivity results, not fit recovery or held-out prediction results.

The stationary task-only example changed the measurement by .12 during [30,60), without changing its underlying latent counterfactual. The numerical stationary x trace varied by roughly 3e-8 around its equilibrium, within the trajectory comparison tolerance. This nuisance ends before the primary persistence window and is not an explanation of that later endpoint. Both full-run figures were visually inspected; labels and legends are readable, and exact trajectory overlap in the confounding figure is expected.

Acceptance is limited to the implemented audit. Optimization-based fit recovery, parameter profiles and complete-schedule held-out model comparison were not implemented; those optional secondary tasks remain explicitly deferred. A mathematical example or an LLM review is not biological validation.
