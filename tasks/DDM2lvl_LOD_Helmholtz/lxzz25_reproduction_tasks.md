Reproduction target: Lu--Xu--Zheng--Zou (2025), Section 5 Tables 5.1--5.9.
Created: 2026-06-01
Updated: 2026-06-03
Verification entry point: `verify/verify_lxzz25_two_level_hybrid_smoke.m`; `verify/verify_lxzz25_article_experiments.m`; `verify/verify_lxzz25_p2_framework.m`; `verify/verify_lxzz25_p2_experiments.m`
Main utilities: `twoLevelHybridSchwarzHelmholtzLOD2D`, `buildLODHelmholtz2D`, `weightedClementP1`, `coarseHatPartition2D`, `assembleHelmholtz2D`, `assemblePlaneWaveBoundaryLoad2D`, `prolongate_P1_P2`, MATLAB `gmres`

# LXZZ25 Two-Level Hybrid Schwarz LOD Reproduction Tasks

This task file records the corrected implementation and verification state for reproducing Section 5 of Lu--Xu--Zheng--Zou, *Two-level hybrid Schwarz Preconditioners for The Helmholtz Equation with high wave number*, arXiv:2408.07669v2.

The corrected instruction file is authoritative:

- `tasks/DDM2lvl_LOD_Helmholtz/codex_instruction_two_level_hybrid_schwarz_corrected_dollar_math.md`

## P2 Fine-Space Non-Paper Extension

The non-paper extension tests the LXZZ hybrid operator on a P2 fine finite element space while keeping the LOD construction in the existing P1 process. The fine mesh size $h$ is the underlying P1 mesh size before P2 extension. For $k \in \{16,32,64,128\}$, the default sweep uses $m=2$, $H=1/k$, and

$$
h^{-1}=\operatorname{align}(\lceil k^{3/2}\rceil,k).
$$

Let $E_{21}$ be the exact P1-to-P2 embedding matrix from `prolongate_P1_P2`. If the existing P1 LOD builder gives $B_1^{\rm trial}$ and $B_1^{\rm test}$, then the P2 coarse basis is

$$
B_2^{\rm trial}=E_{21}B_1^{\rm trial},\qquad
B_2^{\rm test}=E_{21}B_1^{\rm test}.
$$

The P2 coarse matrix is recomputed with the P2 Helmholtz matrix:

$$
A_H=(B_2^{\rm test})^*A_2B_2^{\rm trial},\qquad
A_2=K_2-k^2M_2-ikM_{\Gamma,2}.
$$

The P2 Dirichlet local solver uses P2 DOFs with positive P1 hat weight $w_\ell(x)>0$ as local free DOFs. This removes artificial-boundary DOFs, including artificial/physical-boundary intersections. The P2 impedance local solver uses all local P2 DOFs in the hat support and normalizes partition weights at P2 DOF coordinates.

The P2 extension is verified by `verify/verify_lxzz25_p2_framework.m`. The experiment harness is `verify/verify_lxzz25_p2_experiments.m`, which writes:

- `verify/lxzz25_hybrid_lod/lxzz25_p2_results.csv`
- `verify/lxzz25_hybrid_lod/lxzz25_p2_results.md`

The P2 harness keeps the workstation gates: rows above 300 GB are blocked, rows above 200 GB are marked as requiring permission unless `LXZZ25_P2_ALLOW_GT_200=1`, and default interactive runs are capped by `LXZZ25_P2_MAX_RUN_DOF`.

The crucial correction is that the article's \(Q_m^{(1)}\) and \(Q_m^{(2)}\) are function-level preconditioned operators,
\[
Q_m^{(i)}=B_i^{-1}A,
\]
not residual-to-correction maps. The GMRES left preconditioner must therefore apply \(B_i^{-1}\) to a residual, while `apply` continues to expose \(Q_m^{(i)}\) on coefficient vectors.

## Model and Discretization

The problem is
\[
-\Delta u-\kappa^2u=f\quad\hbox{in }(0,1)^2,\qquad
\partial_nu-i\kappa u=g\quad\hbox{on }\partial\Omega .
\]
The plane-wave exact solution is
\[
u(x,y)=\exp(i\kappa(x+y)/\sqrt2),\qquad
g=i\kappa(d\cdot n-1)u,\quad d=(1/\sqrt2,1/\sqrt2).
\]
For P1 coefficients,
\[
A=K-\kappa^2M-i\kappa M_\Gamma,\qquad D_\kappa=K+\kappa^2M.
\]

The fine mesh follows \(h=O(\kappa^{-3/2})\). With the current default \(C_h=1\), the driver uses
\[
n_h=\lceil C_h\kappa^{3/2}\rceil
\]
and then aligns \(n_h\), \(H\), \(H_{\rm sub}\), and \(\delta\) to integer fine-grid steps.

## Corrected Matrix Form

For the LOD or standard P1 coarse space, the wrapper stores trial and test matrices \(B\) and \(B_\ast\). For LOD,
\[
B=P-C_m,\qquad B_\ast=P-C_m^\ast,
\]
where \(C_m,C_m^\ast\in\mathbb C^{n\times N_H}\) already map coarse coefficients to fine-grid correctors.

The coarse residual inverse and function-level coarse operator are distinct:
\[
M_0^{-1}=B(A_H^{\rm LOD})^{-1}B_\ast^H,\qquad
Q_0=M_0^{-1}A,\qquad
A_H^{\rm LOD}=B_\ast^HAB.
\]
In code:

| Function handle | Meaning |
|---|---|
| `applyM0Inverse(r)` | \(M_0^{-1}r\). |
| `applyQ0(x)` | \(Q_0x=M_0^{-1}Ax\). |
| `applyQ0EuclideanAdjoint(w)` | \(Q_0^H w=A^HB_\ast(A_H^{\rm LOD})^{-H}B^Hw\). |
| `applyEnergyAdjointIMinusQ0(z)` | \((I-Q_0)^{T_D}z=D_\kappa^{-1}(D_\kappa z-Q_0^HD_\kappa z)\). |
| `apply(x)` | Function-level \(Q_m^{(1)}x\) or \(Q_m^{(2)}x\). |
| `applyResidual(r)` | GMRES left preconditioner \(B_1^{-1}r\) or \(B_2^{-1}r\). |

For the Dirichlet variant,
\[
M_D^{-1}=\sum_\ell E_\ell^D A_{\ell,D}^{-1}(E_\ell^D)^H.
\]
For the impedance variant,
\[
M_I^{-1}=\sum_\ell W_\ell C_\ell^{-1}W_\ell^H,\qquad W_\ell=E_\ell X_\ell,
\]
where the local impedance matrix is
\[
C_\ell=K_\ell+\kappa^2M_\ell-i\kappa M_{\partial\Omega_\ell}.
\]

The implemented residual maps are
\[
B_1^{-1}=M_0^{-1}+(I-Q_0)^{T_D}M_D^{-1}(I-AM_0^{-1}),
\]
\[
B_2^{-1}=M_0^{-1}+(I-Q_0)^{T_D}M_I^{-1}(I-AM_0^{-1}).
\]
For a residual \(r\), the corrected update is
\[
z_0=M_0^{-1}r,\qquad r_c=r-Az_0,\qquad z=z_0+(I-Q_0)^{T_D}M_\bullet^{-1}r_c.
\]
The old update \(r-Q_0r\) is incorrect for residual vectors and is no longer used.

## Minimal Verification Matrix

| Check | Status | Verification |
|---|---|---|
| Matrix assembly properties for \(K,M,M_\Gamma,A\) | passed | Smoke Test 1 |
| Plane-wave boundary load consistency under refinement | passed | Smoke Test 2 |
| LOD correctors satisfy \(I_H C_m\approx0\), \(I_H C_m^\ast\approx0\) | passed | Smoke Test 3 |
| LOD patches enlarge monotonically with \(m\) | passed | Smoke Test 4 |
| Correct residual/function identity `applyResidual(A*x) == apply(x)` | passed | Smoke Test 5 |
| Residual update distinguishes \(r-AM_0^{-1}r\) from \(r-Q_0r\) | passed | Smoke Test 5 |
| \((I-Q_0)^{T_D}\) satisfies the \(D_\kappa\)-adjoint identity | passed | Smoke Test 5 |
| Explicit small matrices satisfy \(B_i^{-1}A=Q_m^{(i)}\) | passed | Smoke Test 6 |

Latest smoke output:

```text
========== LXZZ25 Two-Level Hybrid Schwarz Smoke ==========

Test 1: Helmholtz assembly matrix properties ... PASSED
Test 2: Plane-wave consistency under refinement ... PASSED  (5.339e-02 -> 1.437e-02)
Test 3: LOD correctors satisfy Clement kernel constraints ... PASSED  (kernel residual 6.038e-17)
Test 4: LOD patches enlarge monotonically ... PASSED
Test 5: Corrected residual and function-level hybrid identities ... PASSED
Test 6: Explicit small matrices satisfy B_i^{-1} A = Q_m^{(i)} ... PASSED

========== LXZZ25 hybrid smoke tests PASSED ==========
```

## Article Driver and HPC Rule

The article driver writes:

- `verify/lxzz25_hybrid_lod/lxzz25_article_results.csv`
- `verify/lxzz25_hybrid_lod/lxzz25_article_results.md`

Current run settings:

| Setting | Value |
|---|---:|
| Memory limit | 300 GB |
| Per-experiment time limit | 7200 s |
| Interactive run cap | \(N\le 50000\) and coarse elements \(\le 3000\) |
| GMRES tolerance | \(10^{-6}\) |
| GMRES restart | unrestarted |
| Adjoint | \(D_\kappa\)-adjoint |

The 2026-06-01 gated run classified all 76 article rows:

| Status | Count | Meaning |
|---|---:|---|
| `ran` | 4 | Exact smallest rows within the interactive cap reached \(10^{-6}\). |
| `queued_runtime_cap` | 64 | Below 300 GB but above the interactive run cap. |
| `blocked_memory_gt_limit` | 8 | Estimated memory exceeds 300 GB and must not be run. |

Executed exact rows:

| Table | Variant | \(\kappa\) | Paper iterations | Repo iterations | Final relative residual | Status |
|---|---|---:|---:|---:|---:|---|
| 5.1 | \(Q_m^{(1)}\) | 16 | 9 | 15 | \(9.368\times10^{-7}\) | inconsistent |
| 5.2 | \(Q_m^{(2)}\) | 16 | 7 | 8 | \(5.381\times10^{-7}\) | close |
| 5.1 | \(Q_m^{(1)}\) | 32 | 8 | 14 | \(5.591\times10^{-7}\) | inconsistent |
| 5.2 | \(Q_m^{(2)}\) | 32 | 7 | 8 | \(3.956\times10^{-7}\) | close |

The corrected residual-level implementation substantially improves the previous results, especially for \(Q_m^{(2)}\), but the Dirichlet variant is still not paper-matching. The next debugging target is therefore not the residual/function-level matrix split; that now passes explicit tests. The remaining likely causes are the exact article subdomain/partition convention, the localized corrector patch convention, or the authors' precise GMRES algebraic scaling.

## Experiment Parameter Sheet

| Table group | Variant | Parameters | Target |
|---|---|---|---|
| 5.1 | \(Q_m^{(1)}\) | \(h\sim\kappa^{-3/2}\), \(H=\delta=H_{\rm sub}/2\sim\kappa^{-1}\), \(m=\log_2\kappa-1\) | `9,8,8,8,8,8` for \(\kappa=16,32,64,128,256,500\) |
| 5.2 | \(Q_m^{(2)}\) | \(h\sim\kappa^{-3/2}\), \(H\sim\kappa^{-1}\), \(H_{\rm sub}=4H\), \(\delta=2H\), \(m=\log_2\kappa-1\) | `7,7,7,7,7,7` |
| 5.3 | \(Q_m^{(1)}\) | \(\kappa=80\), \(H\sim\kappa^{-1}\), \(m=2\), \(H=\delta=H_{\rm sub}/2\) | `10,9,9,9` |
| 5.4 | \(Q_m^{(2)}\) | \(\kappa=80\), \(H\sim\kappa^{-1}\), \(m=2\), \(H_{\rm sub}=4H\), \(\delta=2H\) | `9,8,8,8` |
| 5.5 | \(Q_m^{(1)}\) | \(\kappa=128\), \(h\sim\kappa^{-3/2}\), \(H=\delta=H_{\rm sub}/2=\kappa^{-1}\) | `8,8,8,8,9,11` |
| 5.6 | \(Q_m^{(2)}\) | \(\kappa=128\), \(h\sim\kappa^{-3/2}\), \(H=\kappa^{-1}\), \(H_{\rm sub}=4H\), \(\delta=2H\) | `7,7,7,7,8,10` |
| 5.7 | \(Q_m^{(1)}\), \(\widetilde Q^{(1)}\) | \(\kappa=40,80,120,160\), \(m=2\), \(H=\delta=H_{\rm sub}/2\) | compare \(H_{\rm sub}=2H_0,H_0\) |
| 5.8 | \(Q_m^{(2)}\), \(\widetilde Q^{(2)}\) | \(\kappa=40,80,120,160\), \(m=2\), \(H=H_0\), \(\delta=H_{\rm sub}/2\) | compare \(\delta=2H_0,4H_0\) |
| 5.9 | \(Q_m^{(1)}\) | fixed fine mesh from \(\kappa_{\max}=120\), \(m=\log_2\kappa\) | small-overlap columns \(\delta=H_0,4h,2h,h\) |

## Open Implementation Notes

The current `coarseHatPartition2D` represents \(\Omega_\ell=\operatorname{supp}\chi_\ell\) using nodal P1 hat supports. This matches the rows where \(\delta=H_{\rm sub}/2\). Table 5.9 varies \(\delta\) while holding \(H_{\rm sub}\) fixed; those rows are currently queued by runtime cap and need a separate small-overlap partition constructor before being treated as reproduced rather than classified.

The exact high-\(\kappa\) rows remain blocked by the 300 GB cap:

- \(\kappa=256,500\) wave-number rows are above 1 TB by the current sparse-LU estimate.
- The fine \(h=2^{-12},2^{-13}\) rows are above 700 GB.
- These must not be launched without explicit permission or a redesigned streamed/distributed local-solve strategy.
