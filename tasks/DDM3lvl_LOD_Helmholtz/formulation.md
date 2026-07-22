Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.
Created: 2026-06-18
Updated: 2026-07-21
Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz.m`; `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`
Main utilities: `twoLevelHybridSchwarzHelmholtz2D`, `buildLODCoarseSchwarzHelmholtz2D`, `buildLODHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`

# Source-Faithful Three-Level LOD-DDM Helmholtz Formulation

## Source-Extracted Hierarchy

Source: `tasks/DDM3lvl_LOD_Helmholtz/Solver_Description.zip`, entry
`main.tex`. The archive contains no PDF file; it contains `main.tex`,
`relerr.png`, and an empty `refs.bib`.

The subscript $\ell$ is a subdomain counter. It is not the overlap or
oversampling parameter. The LOD oversampling parameter is $m$.

The source hierarchy for each coarse-level subdomain is

$$
\operatorname{supp}\chi_{0,\ell}
\subset
\{\chi_{0,\ell}^{>}=1\}
\subset
\operatorname{supp}\chi_{0,\ell}^{>}
\subset
\Omega_{0,\ell},
\qquad
\widetilde\Omega_{0,\ell}
=
\omega_0^{m+C}(\Omega_{0,\ell}).
$$

Here $\Omega_{0,\ell}$ is already an overlapping coarse-level subdomain. It
contains the support of both cutoffs. The local computational domain
$\widetilde\Omega_{0,\ell}$ is not the support of $\chi_{0,\ell}^{>}$; it is
an additional $m+C$ coarse-layer enlargement of $\Omega_{0,\ell}$, where $C$
is a fixed buffer used to keep artificial impedance boundaries away from the
basis functions whose corrected supports are needed inside $\Omega_{0,\ell}$.
The current default is $C=2$ coarse layers.

The safe equality region is the set of coarse nodes/basis functions whose
localized global LOD patches are contained in $\Omega_{0,\ell}$:

$$
S_{0,\ell}^{eq}
=
\left\{
i:
\omega_0^m(\operatorname{supp}\Phi_i)
\subset
\Omega_{0,\ell}
\right\}.
$$

The cutoff construction must satisfy

$$
\{\chi_{0,\ell}^{>}=1\}
\subset
S_{0,\ell}^{eq}.
$$

For columns in this safe region, the locally corrected basis computed on
$\widetilde\Omega_{0,\ell}$ with artificial impedance boundary agrees with the
original localized corrector in the center:

$$
\mathcal C_{\ell,m}\Phi_i
=
\mathcal C_m\Phi_i,
\qquad
i\in S_{0,\ell}^{eq}.
$$

## Local Spaces and Extension

The local fine and coarse spaces are

$$
V_{h,\ell}=V_h(\widetilde\Omega_{0,\ell}),
\qquad
V_{H,\ell}=V_H(\widetilde\Omega_{0,\ell}),
$$

with local quasi-interpolation

$$
I_{H,\ell}:V_h(\widetilde\Omega_{0,\ell})
\to
V_H(\widetilde\Omega_{0,\ell}).
$$

The local fine-scale kernel is

$$
W_{h,\ell}=\ker I_{H,\ell}.
$$

The zero-extension/global coarse-node injection is

$$
E_\ell:V_H(\widetilde\Omega_{0,\ell})\to V_H(\Omega),
\qquad
(E_\ell v)(x)
=
\begin{cases}
v(x), & x\in\widetilde\Omega_{0,\ell}\text{ on the local coarse mesh},\\
0, & \text{otherwise}.
\end{cases}
$$

The local impedance bilinear form on $\widetilde\Omega_{0,\ell}$ is

$$
a_{\varepsilon,\widetilde\ell}(u,v)
=
(\nabla u,\nabla v)_{\widetilde\Omega_{0,\ell}}
-
(k^2+i\varepsilon)(u,v)_{\widetilde\Omega_{0,\ell}}
-
i\sqrt{k^2+i\varepsilon}(u,v)_{\partial\widetilde\Omega_{0,\ell}}.
$$

The local trial and test LOD spaces are

$$
V_{H,m,\ell}
=
(I-\mathcal C_{\ell,m})V_H(\widetilde\Omega_{0,\ell}),
\qquad
V_{H,m,\ell}^{\star}
=
(I-\mathcal C_{\ell,m}^{\star})V_H(\widetilde\Omega_{0,\ell}).
$$

## Source-Mode Local Coarse Solver

The local source-mode matrix is assembled from locally recomputed trial and
test bases on $\widetilde\Omega_{0,\ell}$:

$$
A_{0,\ell}
=
(\Psi_\ell^\star)^*A_{h,\widetilde\ell}\Psi_\ell.
$$

The local coarse solver is the source expression with the local
quasi-interpolation retained:

$$
a_\ell(Q_{0,m,\ell}v_h,v_{H,m,\ell}^{\star})
=
a\!\left(
v_h,
(I-\mathcal C_m^\star)
\Pi_H\bigl(
\chi_{0,\ell}^{>}
E_\ell I_{H,\ell}v_{H,m,\ell}^{\star}
\bigr)
\right).
$$

The gathered additive source-mode coarse preconditioner is

$$
\widetilde M_{0,m}^{-1}r
=
\sum_\ell
\Pi_h\bigl(
\chi_{0,\ell}
E_\ell I_{H,\ell}Q_{0,m,\ell}r
\bigr).
$$

In coefficient form, the current implementation realizes
$E_\ell I_{H,\ell}$ by the local-to-global coarse-node injection and applies
$\chi_{0,\ell}$ and $\chi_{0,\ell}^{>}$ as nodal coarse projection operators.
The partition of unity is required to be exact:

$$
\sum_\ell \Pi_H(\chi_{0,\ell}\cdot)=I
$$

on the coarse coefficient space. This is a hard identity check, not a
tolerance-based floating diagnostic.

This gathered operator is used directly as the default coarse-space solve:

$$
z_0 = \widetilde M_{0,m}^{-1}r_0.
$$

Thus each outer LXZZ preconditioner application uses one application of the
source-DD coarse operator in place of $A_0^{-1}$; it does not run a nested
GMRES solve on $A_0z=r_0$. The preconditioned inner-GMRES mode remains
available only as an explicit diagnostic option. The local source-DD coarse
matrices are stored as sparse LU factors by default; raw sparse backslash
remains available only as an explicit diagnostic option.

## Relation to the LXZZ Wrapper

The fine Helmholtz matrix and positive energy matrix remain

$$
A_h=K_h-M_{k^2+i\varepsilon,h}
-iM_{\sqrt{k^2+i\varepsilon},\Gamma,h},
\qquad
D_h=K_h+M_{|k|^2,h}.
$$

The exact two-level LOD coarse operator remains

$$
A_0=(\Psi^\star)^*A_h\Psi.
$$

The source three-level variant replaces only `coarseSpace.solve` and
`coarseSpace.solveAdjoint` inside `twoLevelHybridSchwarzHelmholtz2D`. The
outer LXZZ residual/function-level algebra and the fine-level impedance local
solver

$$
K_i+k^2M_i-i kM_{\Gamma_i}
$$

are unchanged.

## Current Implementation Stage

Implemented now:

- source-mode coarse geometry with
  $\Omega_{0,\ell}\subset\widetilde\Omega_{0,\ell}$ and
  $\widetilde\Omega_{0,\ell}=\omega_0^{m+C}(\Omega_{0,\ell})$;
- exact nodal coarse partition of unity for $\chi_{0,\ell}$;
- $\chi_{0,\ell}^{>}$ used only on the local test-side restriction and
  enforced to contain $\operatorname{supp}\chi_{0,\ell}$ inside its plateau;
- local source-mode coarse matrices assembled from locally recomputed LOD
  trial/test bases on $\widetilde\Omega_{0,\ell}$ and stored as sparse LU
  factors;
- default one-shot application of the gathered source-DD coarse
  preconditioner in place of $A_0^{-1}$;
- optional parallel construction of independent source-local coarse solvers,
  enabled in the full-sweep runner;
- local-basis comparison classified into safe equality, interior-but-not-safe,
  and artificial-boundary-touched columns;
- injected source coarse apply handles for the existing LXZZ wrapper.

For the iteration-trend sweep, repeated local/global basis comparison is
disabled in the runner after being checked by the focused verification. The
$H_0=1$ one-subdomain control is excluded from the trend table after the
focused verification has checked that it recovers the exact coarse action. This
keeps the $k,H_0,m,\varepsilon$ trend experiment centered on exact two-level
outer GMRES, source three-level outer GMRES using the one-shot coarse DD
operator, and standalone coarse-GMRES diagnostics. Basis fields in that full
CSV are therefore `NaN`.

Postponed to the next diagnostics stage:

- explicit construction of global $M_0^{-1}$;
- explicit construction of $E_0^s$ or $G_0^{(s)}$;
- explicit $|E_0^s|_{G_0}$ and $\alpha_s$ tables for the source geometry.

## Verification Plan

The focused verification script checks:

1. $\operatorname{supp}\chi_{0,\ell}
   \subset \{\chi_{0,\ell}^{>}=1\}\subset\Omega_{0,\ell}$.
2. $\Omega_{0,\ell}\subset\widetilde\Omega_{0,\ell}$ with exactly $m+C$
   coarse expansion layers.
3. $\sum_\ell \Pi_H(\chi_{0,\ell}\cdot)=I$ exactly on coarse coefficients.
4. Safe local/global corrected-basis columns agree to saddle-solve tolerance.
5. `applyM0invAdjoint` satisfies the Euclidean adjoint identity.
6. The LXZZ wrapper still satisfies `applyResidual(A*x) == apply(x)`.
7. Fine-level local solver statistics are unchanged by the coarse-solve
   injection.

The experiment driver records exact two-level GMRES histories, source
three-level GMRES histories, coarse GMRES histories, local-basis errors,
PoU identity status, setup information, local sizes, and memory estimates for
trends in $H_0$ and $k$.
