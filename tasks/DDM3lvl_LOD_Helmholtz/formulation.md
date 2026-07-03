Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.
Created: 2026-06-18
Updated: 2026-06-18
Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz.m`; `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`
Main utilities: `twoLevelHybridSchwarzHelmholtz2D`, `buildLODCoarseSchwarzHelmholtz2D`, `buildLODHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`

# Three-Level LOD-DDM Helmholtz Formulation

## Discrete Objects

The existing LXZZ25 wrapper builds the fine Helmholtz matrix

$$
A_h = K_h - M_{k^2+i\varepsilon,h} - iM_{\sqrt{k^2+i\varepsilon},\Gamma,h},
$$

and the positive energy matrix

$$
D_h = K_h + M_{|k|^2,h}.
$$

For the LOD coarse space, the wrapper stores embedded trial and test bases

$$
\Psi = [\psi_i], \qquad \Psi^\star = [\psi_i^\star].
$$

The exact two-level method uses

$$
A_0 = (\Psi^\star)^*A_h\Psi
$$

and applies the exact coarse residual inverse through `A0 \ r0`. The
three-level method keeps the same $A_0$, $\Psi$, and $\Psi^\star$, but replaces
the exact coarse solve by a one-level Schwarz approximate inverse on the
coarse coefficient space.

The diagnostic norm is never the indefinite Helmholtz matrix. It is

$$
G_0 = \Psi^*D_h\Psi,
\qquad |x|_{G_0}^2=x^*G_0x.
$$

## Coarse Schwarz Approximate Solve

Let $\chi_\ell$ be a coarse partition-of-unity cutoff and let $\chi_\ell^>$
be its enlarged support cutoff. The implementation realizes both as coarse
finite element multiplication matrices, projected back to the nodal coarse
space:

$$
C_\ell=M_H^{-1}M_{\chi_\ell},
\qquad
C_\ell^>=M_H^{-1}M_{\chi_\ell^>}.
$$

For coarse index set $I_\ell$, let $J_\ell$ inject local coarse coefficients
into the global coarse coefficient vector. The local matrix is assembled from
the local fine Helmholtz matrix on the corrected-basis support:

$$
A_{0,\ell}
=
(\Psi_\ell^\star)^* A_{h,\ell}\Psi_\ell.
$$

The approximate inverse is

$$
M_0^{-1}r
=
\sum_\ell
C_\ell J_\ell A_{0,\ell}^{-1}J_\ell^*(C_\ell^>)^*r.
$$

The implemented adjoint action is the exact algebraic adjoint of this
approximate inverse, not a reused primal solve.

## Local-Versus-Global LOD Basis Diagnostic

Before trusting restricted global corrected columns for local coarse matrices,
the smoke test recomputes local LOD trial and test correctors on each coarse
DDM subdomain. For a local subdomain, it builds local fine and coarse meshes,
uses local impedance boundary flags on artificial boundaries, and rebuilds

$$
(I-C_{m,\ell})\Phi_i^\ell,
\qquad
(I-C_{m,\ell}^\star)\Phi_i^\ell.
$$

These local bases are compared with restricted global LOD columns

$$
(I-C_m)\Phi_i,
\qquad
(I-C_m^\star)\Phi_i.
$$

The report separates columns whose global LOD patches are fully contained in
the local domain from columns touching local artificial boundaries. Contained
columns are expected to agree to saddle-solve tolerance. Boundary-touched
columns are diagnostics and can differ because the local problem imposes a
different artificial boundary.

## Diagnostics

The coarse preconditioned and error-propagation operators are

$$
S_0=M_0^{-1}A_0,
\qquad
E_0=I-S_0.
$$

For small coarse systems the implementation assembles `M0inv`, `S0`, and `E0`
explicitly. It then records

$$
|E_0^s|_{G_0}
=
\sqrt{\lambda_{\max}((E_0^s)^*G_0E_0^s,G_0)}
$$

and

$$
\alpha_s
=
\lambda_{\min}\left(
\frac{G_0(I-E_0^s)+(I-E_0^s)^*G_0}{2},
G_0
\right).
$$

The optional $s$-sweep coarse inverse is

$$
G_0^{(s)}r
=
\sum_{j=0}^{s-1}E_0^jM_0^{-1}r,
$$

and satisfies

$$
I-G_0^{(s)}A_0=E_0^s.
$$

## Current Implementation Status

The first implementation is intentionally 2D and P1-only. It preserves the
existing LXZZ25 residual/function-level algebra by injecting a replacement
`coarseSpace.solve` and `coarseSpace.solveAdjoint` into
`twoLevelHybridSchwarzHelmholtz2D`; the fine-level Schwarz solver is not
rewritten.
