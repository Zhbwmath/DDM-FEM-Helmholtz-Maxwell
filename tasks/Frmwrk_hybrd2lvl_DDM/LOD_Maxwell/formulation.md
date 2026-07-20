Reproduction target: LOD coarse Maxwell space inside an LXZZ-type hybrid additive Schwarz solver.
Created: 2026-07-20
Updated: 2026-07-20
Verification entry point: `verify/verify_lod_maxwell_lxzz_hybrid_smoke.m`
Main utilities: `twoLevelHybridSchwarzMaxwell`, `buildMaxwellFineSpace`, `buildLODMaxwellCoarseSpace`, `buildMaxwellDirichletLocalSolvers`, `euclideanTranspose`, `nedelecFreeEdges2D`, `nedelecSubdomainEdges2D`, external reference `LOD4Maxwell/src/LOD/Maxwell/buildLODMaxwell2D`, external reference `LOD4Maxwell/src/LOD/Maxwell/buildLODMaxwell3D`

# LOD Maxwell Coarse Space In LXZZ Hybrid Schwarz

## Plan Changes From The Prior Plan

The integration plan is changed in two ways.

First, the coarse solver uses `LOD4Maxwell` as a reference-only runtime dependency. No `external/` tree, submodule metadata, or copied `LOD4Maxwell` source is added to this repository. The adapter accepts either an already-built `lod` object or an `opts.lodReferencePath` pointing at the sibling checkout. It temporarily appends `LOD4Maxwell/src` to the MATLAB path while building the coarse object, then restores the original path. Because this repo has its own `lodSolveConstrainedSaddle` earlier on the MATLAB path, the adapter defaults reference LOD patch KKT solves to `localSolverMode='direct'` unless the caller explicitly overrides it.

Second, the default adjoint for the Maxwell LXZZ correction is the Euclidean complex adjoint. In MATLAB this is the conjugate transpose operator `'`, recorded in the helper `euclideanTranspose`. Thus the default action of $(I-Q_0)^\dagger$ is $I-Q_0^H$. The energy-adjoint transpose $(I-Q_0)^{T_D}$ is still implemented for comparison, and the wrapper exposes a diagnostic comparing Euclidean, energy, and paper/reference transpose variants.

## Time-Harmonic Maxwell Model

Let $\Omega\subset\mathbb{R}^d$, $d=2$ or $d=3$, and let $\Gamma=\partial\Omega$. The homogeneous Dirichlet, or PEC, time-harmonic Maxwell problem used here is

$$
\nabla\times\nabla\times u-\kappa^2u=f \quad \text{in }\Omega,
$$

with

$$
u\times n=0 \quad \text{on }\Gamma.
$$

In 2D, $u$ is represented by in-plane Nedelec edge fields and $\nabla\times u$ is scalar. In 3D, both $u$ and $\nabla\times u$ are vector fields.

The continuous space is

$$
V=H_0(\operatorname{curl};\Omega)
=\{v\in H(\operatorname{curl};\Omega):v\times n=0 \text{ on }\Gamma\}.
$$

The weak form is: find $u\in V$ such that

$$
a(u,v)=F(v)\qquad \forall v\in V,
$$

where

$$
a(u,v)=(\nabla\times u,\nabla\times v)_{L^2(\Omega)}
-\kappa^2(u,v)_{L^2(\Omega)},\qquad
F(v)=(f,v)_{L^2(\Omega)}.
$$

The algebra below uses coefficient vectors with conjugate-linear testing:

$$
a(u_h,v_h)=v^H A_h u,\qquad F(v_h)=v^H b_h.
$$

## Fine Edge Space

Let $\mathcal{T}_h$ be the fine triangular or tetrahedral mesh and let $E_h$ be its oriented edge set. The lowest-order Nedelec space is

$$
V_h^{\rm full}=\mathcal{N}_1(\mathcal{T}_h).
$$

Boundary edge DOFs are

$$
E_h^\Gamma=\{e\in E_h:e\subset\Gamma\}.
$$

The active Dirichlet edge DOF set is

$$
E_h^0=E_h\setminus E_h^\Gamma,
$$

and the algebraic fine space is

$$
V_h=\{v_h\in V_h^{\rm full}: \operatorname{dof}_e(v_h)=0\ \forall e\in E_h^\Gamma\}
\cong \mathbb{C}^{N_h},\qquad N_h=|E_h^0|.
$$

Let $S_h^{\rm full}$ be the curl-curl matrix and $M_h^{\rm full}$ the Nedelec mass matrix on all edges:

$$
(S_h^{\rm full})_{ij}=(\nabla\times\phi_j,\nabla\times\phi_i),\qquad
(M_h^{\rm full})_{ij}=(\phi_j,\phi_i).
$$

Let $R_h$ select the active edge rows $E_h^0$. The reduced fine matrices are

$$
S_h=R_hS_h^{\rm full}R_h^T,\qquad
M_h=R_hM_h^{\rm full}R_h^T,
$$

and

$$
A_h=S_h-\kappa^2M_h,\qquad
D_h=S_h+\kappa^2M_h.
$$

$D_h$ is the positive energy matrix used only to define the optional energy adjoint and diagnostics.

## LOD Coarse Space

Let $\mathcal{T}_H$ be a nested coarse mesh. The coarse active edge space is

$$
V_H=\mathcal{N}_1(\mathcal{T}_H)\cap H_0(\operatorname{curl};\Omega)
\cong\mathbb{C}^{N_H}.
$$

The nested Nedelec interpolation is

$$
P_{\operatorname{curl}}:\mathbb{C}^{N_H}\rightarrow\mathbb{C}^{N_h}.
$$

Let $S_h^{\rm sc}$ be the scalar P1 Dirichlet stiffness matrix, $M_h^{\rm sc}$ the scalar P1 Dirichlet mass matrix, and

$$
G_h:\mathbb{C}^{n_h^{\rm sc}}\rightarrow\mathbb{C}^{N_h}
$$

the discrete gradient from scalar nodal P1 DOFs to Nedelec edge DOFs. Let

$$
R_{\rm sc}:\mathbb{C}^{n_H^{\rm sc}}\rightarrow\mathbb{C}^{n_h^{\rm sc}}
$$

be the nested scalar P1 interpolation. The coarse scalar mass-moment constraint matrix is

$$
C_{\rm sc}=M_h^{\rm sc}R_{\rm sc}.
$$

The scalar-gradient corrector potential $Q_{\rm grad}$ is defined by the constrained saddle problem

$$
\begin{bmatrix}
S_h^{\rm sc} & C_{\rm sc}\\
C_{\rm sc}^H & 0
\end{bmatrix}
\begin{bmatrix}
Q_{\rm grad}\\ \Lambda_{\rm sc}
\end{bmatrix}
=
\begin{bmatrix}
G_h^HM_hP_{\operatorname{curl}}\\ 0
\end{bmatrix}.
$$

The lifted scalar-gradient corrector is

$$
K_{\rm grad}=G_hQ_{\rm grad}.
$$

Let $\mathcal{I}_{\rm nc}$ denote fine scalar Dirichlet nodes that are not coarse scalar nodes. The constraint span for the $X$-corrector is

$$
X_r=\operatorname{range}
\begin{bmatrix}
P_{\operatorname{curl}} & G_h(:,\mathcal{I}_{\rm nc})
\end{bmatrix},
$$

and its mass constraint matrix is

$$
J_X=M_h
\begin{bmatrix}
P_{\operatorname{curl}} & G_h(:,\mathcal{I}_{\rm nc})
\end{bmatrix}.
$$

The $X$-corrector $W_X$ is defined by

$$
\begin{bmatrix}
A_h & J_X\\
J_X^H & 0
\end{bmatrix}
\begin{bmatrix}
W_X\\ \Lambda_X
\end{bmatrix}
=
\begin{bmatrix}
A_hP_{\operatorname{curl}}\\ 0
\end{bmatrix}.
$$

The corrected Maxwell LOD coarse injection is

$$
P_0=P_{\operatorname{curl}}-W_X-G_hQ_{\rm grad}.
$$

For the reference-only implementation, $P_0$ is reconstructed from the `LOD4Maxwell` output fields instead of importing the source files.

The coarse matrix and solve are

$$
A_0=P_0^HA_hP_0,\qquad
M_0^{-1}=P_0A_0^{-1}P_0^H.
$$

This first implementation uses the same trial and test corrected space, so $P_0^\ast=P_0$.

## Local Dirichlet Schwarz Solver

Let $\{\Omega_\ell\}_{\ell=1}^{N_s}$ be overlapping subdomains from `partitionMesh2D` or `partitionMesh3D`. Each local solve uses homogeneous Dirichlet conditions on the artificial local boundary.

For each subdomain, define the local interior edge set

$$
E_\ell^0
=\{e\in E_h^0:\operatorname{star}(e)\subset\Omega_\ell\},
$$

where $\operatorname{star}(e)$ is the set of fine elements incident to edge $e$. Let

$$
R_\ell:\mathbb{C}^{N_h}\rightarrow\mathbb{C}^{|E_\ell^0|}
$$

select the reduced fine DOFs belonging to $E_\ell^0$, and let

$$
E_\ell=R_\ell^H
$$

be the zero-extension. The local Dirichlet Maxwell matrix is

$$
A_\ell=R_\ell A_hE_\ell.
$$

The one-level local additive Schwarz inverse is

$$
M_{\rm loc}^{-1}=\sum_{\ell=1}^{N_s}E_\ell A_\ell^{-1}R_\ell.
$$

The implementation stores sparse LU factors by default because $A_h$ is indefinite for time-harmonic Maxwell.

## LXZZ Hybrid Operator

The coarse projection is

$$
Q_0=M_0^{-1}A_h=P_0A_0^{-1}P_0^HA_h.
$$

The function-level LXZZ two-sided hybrid operator is

$$
Q_m=Q_0+(I-Q_0)^\sharp M_{\rm loc}^{-1}A_h(I-Q_0),
$$

where $\sharp$ is the selected transpose or adjoint convention.

The residual preconditioner used by GMRES is

$$
B^{-1}=M_0^{-1}+(I-Q_0)^\sharp M_{\rm loc}^{-1}(I-A_hM_0^{-1}).
$$

The verifier checks the operator identity

$$
B^{-1}A_hx=Q_mx
$$

on deterministic complex vectors.

## Transpose And Adjoint Choices

The default Maxwell choice is the Euclidean complex adjoint:

$$
Q_0^H=A_h^HP_0A_0^{-H}P_0^H.
$$

Therefore

$$
(I-Q_0)^\sharp=(I-Q_0)^H=I-Q_0^H.
$$

In MATLAB this uses `euclideanTranspose(A)`, which returns `A'`, the complex-conjugate transpose. It deliberately does not use `A.'`.

The optional energy adjoint is

$$
Q_0^{T_D}=D_h^{-1}Q_0^HD_h,
\qquad
(I-Q_0)^{T_D}=D_h^{-1}(I-Q_0)^HD_h.
$$

The paper/reference transpose kept for comparison is

$$
\widetilde Q_0^H v=P_0A_0^{-H}P_0^HA_h^Hv.
$$

It is not the Euclidean matrix adjoint of $Q_0$ unless the operators commute in a special way. The smoke verifier reports the relative differences between the Euclidean default and the energy/paper variants.

## Implementation Plan

1. Add 2D Nedelec boundary, free-edge, edge-incidence, and subdomain-edge utilities mirroring the existing 3D utilities.
2. Add `buildMaxwellFineSpace` for the reduced Dirichlet matrix $A_h=S_h-\kappa^2M_h$ and energy matrix $D_h=S_h+\kappa^2M_h$ in 2D and 3D.
3. Add `buildLODMaxwellCoarseSpace` as a reference-only adapter for the sibling `LOD4Maxwell` checkout, with an already-built `lod` object as an alternative input.
4. Add `buildMaxwellDirichletLocalSolvers` for local edge-star Dirichlet subdomain solves.
5. Add `twoLevelHybridSchwarzMaxwell` implementing the LXZZ function-level and residual-level actions with `adjointType='euclidean'` as default.
6. Add a focused smoke verifier that covers 2D and 3D reduced spaces, LOD coarse basis normalization, local Dirichlet edge sets, Euclidean versus energy/paper transpose diagnostics, and GMRES application.

## Scope Limits

This phase supports constant scalar $\kappa$ and homogeneous Dirichlet/PEC boundary conditions. Variable coefficients, impedance local solvers, and direct inclusion of `LOD4Maxwell` files in this repository are out of scope.
