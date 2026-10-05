# Analytical constraints and interpretation

Evidence class: exact properties of specified dimensionless equations. These are not empirical findings. Prepared before the v1 full simulation run, 5 October 2026.

## Bounds and fixed systems

For M0, M1 and M2, the vector field at x=0 is a>=0 and at x=1 is -r<=0. For M4 the latter is -r(1-d)<=0. At d=0, its derivative is nonnegative; at d=1 it is -eta<=0. Nonnegative commanded drive keeps the M3 input state nonnegative. Locally Lipschitz dynamics therefore preserve the admissible state domain. This is a mathematical property; numerical implementations must still check rather than clip violations.

For a,r>0, M0 has the unique attracting equilibrium x*=a/(a+r). M3 has the same settled x* with s*=a; its Jacobian is triangular with eigenvalues -(a+r) and -1/tau. Different finite-time endpoints or a ramp loop in M3 do not demonstrate two attractors.

For M2 with a>0 and b>0, the quadratic equilibrium polynomial is -b*x^2+(b-a-r)*x+a. Its root product is -a/b<0, so exactly one root is positive. With r>0, the signs at 0 and 1 place that root strictly inside the state interval. The derivative there equals minus the positive square root of the discriminant, proving attraction. Degenerate zero-rate cases require separate boundary handling.

## M1 folds and history

Write f(x)=a(1-x)+b*x^2*(1-x)-r*x. At a fold, f=0 and f_x=0 imply

`a=b*x^2*(1-2*x)` and `r=2*b*x*(1-x)^2`.

For b>0 and nonnegative a, the relevant fold parametrization has 0<=x<=1/2. The cusp occurs at x=1/3, a=b/27 and r=8*b/27. Away from folds, three physical simple roots with alternating negative, positive, negative derivatives establish two attracting equilibria and their separating unstable root. A numerical near-fold root match alone is insufficient to certify settling.

At (a,b,r)=(.06,4,1.05), the roots are (3-sqrt(5.25))/10, .4 and (3+sqrt(5.25))/10, approximately .070871215, .4 and .529128785. Under the pretreatment drive a=.14 the low attractor is absent. A successful low-state pulse under maintained a=.06 therefore need not persist after all controls restore a=.14. The two policies describe different final dynamical systems.

## Imposed treatment ordering

For scalar M0-M2 and the scalar M1 frozen-d variant, increasing r changes the derivative by -delta_r*x<=0; decreasing a changes it by delta_a*(1-x)<=0. From a common initial state, the scalar comparison principle orders trajectories under pointwise ordered controls. Under identical reduced drive, an added nonnegative repair pulse cannot increase x. The sign of this benefit is built into the equation. Its magnitude, target crossing, dependence on final policy and persistence are informative conditional results. A factorial contrast on an arbitrary x scale does not establish pharmacological synergy.

For M3, s is independent of repair and order-preserving under ordered input commands; the corresponding x comparison follows. Extending scalar reasoning to evolving M4 requires a separate coupled-system argument and must not be asserted solely by analogy.

## Structural and observation restrictions

If eta=0, d is nondecreasing. An initial d above a specified structural target cannot meet that target later within this model. If eta=epsilon=0, d is a fixed parameter. When eta=0 and x<=xc, evolving M4 can have neutral directions in d; frozen-d x stability does not make the full equilibrium an isolated attractor. The proposed readout q=(1-x)(1-d) satisfies q<=1-d; it is neither calibrated cognition nor a patient score.

With frozen unknown d and no independently known additive repair perturbation, x dynamics depend on r and d only through rho=r(1-d). Pairs sharing rho are exactly indistinguishable from x alone under the same a,b and drive schedule. Known additive changes to r can break this particular algebraic symmetry; confounding must be stated for the actual control/measurement design.

At equilibrium, multiplying all of a,b,r by the same positive factor preserves roots while changing time scales. Steady-state data alone cannot identify that common rate scale. Dynamic observations with an anchored clock can contain additional information.

The instantaneous state Jacobians of y1=1-x and y2=.85(1-x) are collinear. Adding an anchored structural readout y2=d increases their instantaneous rank. This rank statement is not a complete dynamic observability or unknown-parameter identifiability proof. Identical raw nuisance shifts in two tasks with unequal known loadings can be distinguishable from x; a nuisance proportional to those loadings is confounded. The observation experiment must keep these cases separate.

## Recovery and cost claims

A successful tested regimen is a feasibility witness for that regimen family, state, target and horizon. Failure of a finite amplitude-duration grid does not prove universal inability to recover. Normalized quadratic control expenditure is arbitrary exposure accounting, not toxicity, money or clinically optimal dose. A grid frontier is a frontier among sampled controls, not a solved continuous optimal-control problem.
