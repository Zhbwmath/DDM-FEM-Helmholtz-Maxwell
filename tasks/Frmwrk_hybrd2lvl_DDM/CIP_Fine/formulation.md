Reproduction target: CIP fine-space instance of the abstract two-sided two-level hybrid DDM.
Created: 2026-06-26
Updated: 2026-07-03
Verification entry point: `verify/verify_cip_lxzz_lod_medium.m`
Main utilities: `assembleHelmholtzCIP2D`, `assembleCIP2D`, `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildLODHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`

# CIP Fine-Space Instance

## Model Problem

The current two-dimensional scalar Helmholtz instance uses $\Omega=(0,1)^2$ and wave number $k>0$:

$$
-\Delta u-k^2u=f \quad \text{in }\Omega .
$$

The impedance boundary convention used by the assembled fine matrix is

$$
\partial_n u - i k u = g \quad \text{on }\partial\Omega ,
$$

so the boundary contribution in the homogeneous operator is $-ik\int_{\partial\Omega}u\bar v\,ds$.

## Fine Space

For the default medium verification, $V_h$ is the continuous P1 Lagrange space on a uniform triangular mesh of size $h$. The same framework path now supports P2 and P3 fine spaces through the existing high-order Lagrange assemblers and P1-to-P2/P1-to-P3 injection matrices. The wave-resolution rule is

$$
h^{-1}=\operatorname{align}\left(\left\lceil C_h k^{(2p+1)/(2p)}\right\rceil,k\right),
\qquad p=1,\qquad C_h=1,
$$

where `align` rounds upward so the fine mesh is nested with the coarse scale $H=1/k$.

## CIP Sesquilinear Form

The fine-level sesquilinear form is

$$
a_{\mathrm{CIP}}(u,v)
=(\nabla u,\nabla v)_\Omega
-k^2(u,v)_\Omega
-ik(u,v)_{\partial\Omega}
+J_\gamma(u,v).
$$

For polynomial degree $p$, the continuous interior penalty term is

$$
J_\gamma(u,v)
=\sum_{e\in\mathcal{E}_h^{\mathrm{int}}}
\sum_{j=1}^{p}
i\gamma_j h_e^{2j-1}
\int_e [\partial_n^j u]\,
\overline{[\partial_n^j v]}\,ds .
$$

In the current P1 run only the $j=1$ normal-derivative jump is active. The project assembler supports the same contract for higher Lagrange degrees, but this subtask starts with the P1 case.

With coefficient vectors,

$$
A_{\mathrm{CIP}}=K-k^2M-ikM_\Gamma+J_\gamma .
$$

The $D_h$ energy inner product for the hybrid adjoint is the positive Helmholtz energy matrix

$$
D_h=K+k^2M .
$$

The penalty contribution $J_\gamma$ remains in the fine operator $A_{\mathrm{CIP}}$ and local matrices, not in $D_h$.

## Coarse Space

This instance deliberately uses the normal-FEM LOD/Helmholtz-harmonic coarse space rather than a CIP-generated coarse basis. Let $\Psi$ and $\Psi^\ast$ denote the normal-FEM LOD trial and test basis matrices produced from the auxiliary Helmholtz FEM operator on the coarse mesh.

The coarse injections for the abstract framework are

$$
P_0=E_h\Psi,\qquad P_0^\ast=E_h\Psi^\ast ,
$$

where $E_h$ injects the normal P1 nodal basis into the active fine space. For P1 CIP, $E_h$ is the identity on nodal DOFs. For LOD P1-P2, the coarse space and correctors remain P1, then the corrected P1 basis is injected into the P2 fine space before the coarse matrix is recomputed. The coarse matrix is always assembled against the active CIP fine operator:

$$
A_0=(P_0^\ast)^H A_{\mathrm{CIP}} P_0 .
$$

This is the intended separation: normal-FEM coarse basis construction, CIP fine-level Galerkin projection.

## Local Solvers

Both local solver variants are tested under the same two-sided hybrid algebra.

The Dirichlet variant uses overlapping or hat-support subdomain restrictions with local matrices extracted from the CIP Helmholtz form on the local free DOFs:

$$
M_{\mathrm{loc,D}}^{-1}
=\sum_j E_j A_{\mathrm{CIP},j}^{-1}R_j .
$$

The impedance variant uses local CIP Helmholtz matrices with Robin terms on artificial subdomain interfaces:

$$
a_{\mathrm{loc,I},j}(u,v)
=a_{\mathrm{CIP},j}(u,v)
-ik(u,v)_{\Gamma_j^{\mathrm{art}}}.
$$

The resulting local action has the same abstract form

$$
M_{\mathrm{loc,I}}^{-1}
=\sum_j E_j A_{\mathrm{loc,I},j}^{-1}R_j .
$$

The extension $E_j$ includes the project partition-of-unity or hat weighting used by the current local-solver builder.

## Two-Sided Hybrid Preconditioner

For either local variant,

$$
Q_0=M_0^{-1}A_{\mathrm{CIP}},\qquad
M_0^{-1}=P_0A_0^{-1}(P_0^\ast)^H ,
$$

and the residual GMRES preconditioner is

$$
B^{-1}
=M_0^{-1}
+(I-Q_0)^{T_D}M_{\mathrm{loc}}^{-1}(I-A_{\mathrm{CIP}}M_0^{-1}).
$$

The corresponding function-level operator is

$$
Q_m
=Q_0+(I-Q_0)^{T_D}M_{\mathrm{loc}}^{-1}A_{\mathrm{CIP}}(I-Q_0).
$$

The verifier checks $B^{-1}A_{\mathrm{CIP}}=Q_m$ before running GMRES.

## Medium Sweep Parameters

The first sweep uses:

- $k\in\{16,32,64,128\}$;
- $p=1$ and $C_h=1$;
- $H=1/k$;
- LOD oversampling $m=\max(1,\operatorname{round}(\log_2 k-1))$ for $k<128$, and fixed $m=3$ for $k\ge 128$ to reduce high-frequency memory cost;
- local variants `dirichlet` and `impedance`, with LXZZ local partition spacing $1/k$ for Dirichlet and $2/k$ for impedance;
- local solver mode `lu`, LOD solver mode `direct`, energy-adjoint hybrid, and parfor-enabled subdomain setup as in the LXZZ article driver;
- local storage `matrix`, meaning local sparse matrices are stored at setup time and factored inside each GMRES preconditioner application instead of storing all local LU factors; local apply mode defaults to `auto`, which keeps the faster full-vector worker outputs below the configured memory threshold and switches to compact local contributions above it; by default, high-$k$ local setup and LOD construction are serial to avoid parfor worker serialization failures, while local solves are parallelized during the preconditioner application;
- GMRES tolerance $10^{-6}$ and maximum iteration count 100.

The force-run command for the high-$k$ rows is:

```powershell
$env:CIP_LXZZ_LOD_KVALUES='64,128'
$env:CIP_LXZZ_LOD_FORCE_RUN='1'
$env:CIP_LXZZ_LOD_PARFOR='1'
$env:CIP_LXZZ_LOD_PARPOOL_WORKERS='8'
$env:CIP_LXZZ_LOD_LOCAL_SETUP_PARFOR='0'
$env:CIP_LXZZ_LOD_LOCAL_APPLY_PARFOR='1'
$env:CIP_LXZZ_LOD_LOD_PARFOR='0'
$env:CIP_LXZZ_LOD_HIGHK_OVERSAMPLING='3'
$env:CIP_LXZZ_LOD_LOCAL_STORAGE='matrix'
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_cip_lxzz_lod_medium.m');"
```
