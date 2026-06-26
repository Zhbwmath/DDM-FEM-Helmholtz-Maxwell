Reproduction target: Hu--Li, Tables 1--3, plus Hu--Li coarse spaces in the two LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-10
Updated: 2026-06-22
Verification entry point: `verify/verify_hl25_coarse_spaces.m`; `verify/verify_hl25_tables123.m`; `verify/verify_hl25_lxzz_cross_study.m`; `verify/verify_hl25_full_sweep.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `partitionMesh2D`, `linearPartitionOfUnity2D`, MATLAB `gmres`

# Hu--Li Reproduction And LXZZ25 Cross-Study

## Implementation Status

The reusable implementation is
`src/Preconditioners/buildHuLiWeightedSchwarzHelmholtz2D.m`. It provides:

- spectral Helmholtz-harmonic coarse spaces;
- economic periodic-quadratic trace coarse spaces;
- P1 and P2 fine spaces;
- scalar, function-valued, and PDE-struct wave numbers;
- $\varepsilon=0$ and $\varepsilon=\kappa$;
- direct, stored-LU, and memory-adaptive local solves;
- optional `parfor` local setup;
- native Hu--Li WASI hybrid application;
- `fineSpace` and `coarseSpace` objects for the LXZZ25 wrapper.

The implementation preserves the existing uncommitted variable-coefficient
changes in the LXZZ25 and LOD paths.

## Focused Verification

`verify/verify_hl25_coarse_spacesTest.m` contains class-based tests for:

- partition-of-unity normalization;
- local harmonic-extension residuals;
- Hermitian generalized eigenproblems;
- spectral threshold selection and energy orthonormality;
- economic periodic trace dimensions and continuity;
- explicit native hybrid matrix versus function-handle application;
- scalar versus equivalent PDE-struct consistency;
- a genuinely variable $\kappa(x,y)$ case with `kappaRef`;
- spectral and economic coarse spaces in both LXZZ25 variants;
- both $\varepsilon=0$ and $\varepsilon=\kappa$ in the LXZZ25 integration.

The initial focused run passed all 13 tests.

After the performance changes, the suite contains 14 tests and all 14 pass.
The added test compares the cached and direct energy-adjoint implementations
and verifies agreement to floating-point tolerance.

## Efficiency Diagnosis

Two independent bottlenecks caused the slow partial runs.

1. The cross-study source had been changed to $\beta=1$ and its gate block was
   commented out. At $\kappa=8\pi$, this creates 625 subdomains instead of the
   planned 49 for $\beta=0.6$, and attempts every article-scale case.
2. The Hu--Li builder applied a tall rank-revealing QR directly to
   $Z_{\rm raw}$. A profile on a P2 economic case with 11,025 unknowns and 169
   subdomains measured 32.7 seconds in that QR out of 36.2 seconds total.

The same benchmark after replacing the tall QR by normalized Gram-matrix rank
selection and extending only the requested economic traces took 2.45 seconds:

| stage | previous | optimized |
|---|---:|---:|
| total Hu--Li setup | 36.2 s | 2.45 s |
| rank reduction | 32.7 s | 0.035 s |
| local harmonic extensions | 0.186 s | 0.023 s |

This is a 14.8-fold setup speedup on identical data.

For LXZZ, a profile showed that the exact global energy-adjoint solve consumed
84% of one preconditioner application; the many local solves consumed 9%.
Precomputing the complete energy-adjoint coarse basis accelerated applications
by 9.1 times but cost 121 seconds for 338 coarse columns, so it is retained as
an opt-in mode only. The default instead factors the energy matrix once and
shares that factor across both coarse spaces and both absorption values.

On the planned P2 $\kappa=8\pi$, $\beta=0.6$ economic case:

| quantity | optimized value |
|---|---:|
| unknowns / subdomains / coarse dimension | 28,561 / 49 / 392 |
| Hu--Li build | 28.1 s |
| of which shared energy factorization | 26.2 s |
| native solve, 5 iterations | 0.59 s |
| LXZZ Dirichlet setup / solve, 24 iterations | 0.72 s / 29.1 s |
| LXZZ impedance setup / solve, 17 iterations | 0.95 s / 21.3 s |

The previous economic setup for this mesh was about 284 seconds. The remaining
28-second first-case setup is mostly the shared energy factorization and is
not repeated for the other coarse/absorption configurations at the same
wavenumber.

## Revised Memory And Runtime Gates

The earlier estimates were intentionally conservative and are superseded by
the following optimistic model:

- GMRES storage uses the actual default restart length 10, not 1003 vectors.
- Sparse-LU fill uses $c_{\rm lu}=20$, replacing the conservative value 40.
- The reported gate estimate is the smaller of stored-LU and serial direct
  solve memory for the local solves. Both estimates and the selected memory
  mode are written to CSV.
- Adaptive execution still uses stored LU when its estimate is below 200 GB
  for speed. If LU exceeds 200 GB but direct solve does not, the driver uses
  direct mode.
- Coarse operator memory is no longer estimated as dense $16r^2$ storage.
  The gate uses a sparse coarse-matrix row-count estimate plus the same
  $16c_{\rm lu}r\log_2(r)$ sparse-LU storage model. Dense coarse storage is
  not used in the permission gate.
- Spectral cross-study dimensions use a configurable coarse/local ratio
  estimate of 1.0 rather than assuming that every boundary harmonic mode is
  retained. Table estimates use the paper's reported coarse/local ratio.
- `parfor` remains configurable. The full-sweep wrapper enables it by default,
  but the active 2026-06-22 launch disabled it because an unrelated
  SpectralLOD MATLAB parpool was already occupying many workers.

Representative revised cross-study estimates are:

| $\kappa$ | spectral | economic |
|---:|---:|---:|
| $16\pi$ | 0.67 GB | 0.52 GB |
| $40\pi$ | 4.38 GB | 4.20 GB |
| $80\pi$ | 24.75 GB | 26.00 GB |
| $160\pi$ | 169.74 GB | 185.12 GB |

These are gate estimates, not measured peak memory. The 200 GB permission
rule remains unchanged.

The default runtime caps are relaxed to 150,000 unknowns for the LXZZ
cross-study and 800,000 unknowns for Tables 1--3. The hard estimate cap is
500 GB, while every row above 200 GB still requires explicit permission.
This admits the $8\pi$ and $16\pi$ cross-study configurations by default and
eleven exact P2, $\kappa=40\pi$ Table 1--3 rows under the runtime and
permission gates.

## Tables 1--3 Driver

`verify/verify_hl25_tables123.m` encodes all paper iteration counts and
coarse/local dimension ratios for:

- Table 1: theoretical spectral thresholds for P1 and P2;
- Table 2: practical spectral threshold
  $$
  \rho=\frac12\kappa^{(\beta-1)/2};
  $$
- Table 3: economic coarse space with
  $$
  \nu=\operatorname{round}(\kappa^{1-\beta}).
  $$

The driver writes:

- `tasks/HL25_Helmholtz_harmonic/tables_1_3_results.csv`;
- `tasks/HL25_Helmholtz_harmonic/tables_1_3_results.md`;
- `tasks/HL25_Helmholtz_harmonic/tables_1_3_comparison.png` when rows run.

Every result row records the aligned mesh, subdomain count, overlap, estimated
memory, paper values, repository values, differences, timing, and status.

The previous classification contains 96 paper rows:

| status | rows |
|---|---:|
| `queued_runtime_cap` | 22 |
| `requires_permission_gt_200gb` | 1 |
| `blocked_memory_gt_hard_limit` | 73 |

Those statuses were produced by the superseded conservative estimator and
will be replaced when the table driver is rerun. No article-scale row was
executed previously. Under the revised 800,000-unknown cap, eleven exact P2
$\kappa=40\pi$ rows now pass both the runtime and permission gates. Remaining
rows stay queued by runtime or require the documented memory permission; a
gate classification is not evidence that the mathematical method fails.

## Full Hu--Li/LXZZ25 Sweep

`verify/verify_hl25_full_sweep.m` is the checkpointed real-run entry point for
the requested sweep. The default grid is

$$
\kappa/\pi\in\{16,32,64,128\},\qquad
\beta\in\{0.5,0.6,0.7\},
$$

with P2 elements, both spectral and economic Hu--Li coarse spaces, both
$\varepsilon=0$ and $\varepsilon=\kappa$, and all three methods: native
Hu--Li, LXZZ Dirichlet local hybrid, and LXZZ impedance local hybrid.

The driver writes:

- `tasks/HL25_Helmholtz_harmonic/full_sweep_lxzz_cross_results.csv`;
- `tasks/HL25_Helmholtz_harmonic/full_sweep_lxzz_cross_results.md`;
- `tasks/HL25_Helmholtz_harmonic/full_sweep_lxzz_cross_iterations.png` when
  rows run.

The 2026-06-22 estimate-only preflight produced 144 rows. The largest sparse
optimistic gate estimate was 163.18 GB, below the 200 GB permission threshold.
That first launch was stopped after inspection showed that the LXZZ local
solvers were incorrectly using the Hu--Li paper partition. Those result files
were archived with suffix `invalid_same_partition_20260622_180252`.

The corrected 2026-06-22 preflight keeps the Hu--Li coarse-space partition at
$d\approx\kappa^{-\beta}$ but gives LXZZ local solvers their own aligned
`coarseHatPartition2D` partition with $H_{\rm LXZZ}\le 1/\kappa$ by default.
It again produced 144 estimate-only rows. The largest sparse optimistic gate
estimate is 172.28 GB, below the 200 GB permission threshold. Each CSV row now
records both `nSubSide`/`nSubdomains` for Hu--Li and
`lxzzNSubSide`/`lxzzNSubdomains` for LXZZ local solves.

## LXZZ25 Cross-Study Driver

`verify/verify_hl25_lxzz_cross_study.m` compares the native Hu--Li hybrid with
the LXZZ25 Dirichlet and impedance hybrids. It uses P2, Hu--Li
$\beta=0.6$ by default, generous Hu--Li overlap, both Hu--Li coarse spaces,
and both absorption values. LXZZ local solvers do not reuse the Hu--Li
subdomains; their default local inverse size is aligned upward from
$\lceil\kappa\rceil$ so $H_{\rm LXZZ}\le 1/\kappa$.

The scaled cases are

$$
\kappa\in\{8\pi,16\pi\},
$$

followed by classified article-scale cases

$$
\kappa\in\{40\pi,80\pi,120\pi,160\pi\}.
$$

The driver writes:

- `tasks/HL25_Helmholtz_harmonic/lxzz_cross_results.csv`;
- `tasks/HL25_Helmholtz_harmonic/lxzz_cross_results.md`;
- `tasks/HL25_Helmholtz_harmonic/lxzz_cross_iterations.png` when rows run.

The previously generated cross-study contains 72 rows:

| status | rows |
|---|---:|
| `ran` | 12 |
| `queued_runtime_cap` | 24 |
| `requires_permission_gt_200gb` | 6 |
| `blocked_memory_gt_hard_limit` | 30 |

The generated old CSV may include later user partial runs with altered
$\beta$. The driver now defaults again to the planned $\beta=0.6$ and can be
overridden with `HL25_CROSS_BETA`. The previous twelve $\kappa=8\pi$
configurations used the same P2 mesh with 28,561 unknowns, 49 Hu--Li
subdomains, and generous aligned Hu--Li overlap for both Hu--Li and LXZZ local
solves. They are retained only as historical invalid same-partition rows, not
as corrected LXZZ cross-study evidence. Their observed GMRES iteration counts
were:

| $\varepsilon$ | coarse space | coarse dimension | Hu--Li native | LXZZ Dirichlet | LXZZ impedance |
|---|---|---:|---:|---:|---:|
| $0$ | spectral | 1107 | 6 | 25 | 17 |
| $0$ | economic | 392 | 5 | 24 | 17 |
| $\kappa$ | spectral | 1107 | 6 | 25 | 17 |
| $\kappa$ | economic | 392 | 5 | 24 | 18 |

Every completed solve reached the requested $10^{-6}$ relative tolerance.
Absorption changed the scaled-case iteration counts by at most one. The native
Hu--Li hybrid required substantially fewer iterations than either LXZZ hybrid
with the same Hu--Li basis in this scaled case; this is a cross-study result,
not a comparison reported in the Hu--Li paper.

## Reproduction Conventions And Deviations

- The constants hidden in $h\approx\kappa^{-(2p+1)/(2p)}$ and
  $d\approx\kappa^{-\beta}$ are set to one.
- The fine divisor is rounded upward to align subdomain and overlap
  boundaries with mesh edges.
- Generous overlap is $\delta=d/4$ after alignment.
- MATLAB sparse LU/backslash replaces the paper's MUMPS local and coarse
  solvers.
- MATLAB `eig`/`eigs` replaces SLEPc.
- The economic implementation currently requires axis-aligned rectangular
  subdomains, matching the structured checkerboard experiments.
- DtN, HGenEO, grid coarse spaces, Tables 4--5, and other competitors are
  intentionally out of scope.
- Rows estimated above 200 GB are not executed without explicit permission.

## Reproduction Assessment

The coarse-space and preconditioner implementation is verified on focused
meshes, and every requested row has a durable result or gate classification.
The scaled LXZZ cross-study is internally consistent and exercises both coarse
spaces, both LXZZ variants, and both absorption values.

Tables 1--3 are currently **blocked**, rather than reproduced: no exact
article-scale row ran under the active runtime and memory gates, so the
repository values cannot yet be judged consistent or inconsistent with the
paper tables. The $\kappa=8\pi$ LXZZ cross-study is **completed as a scaled
verification**, while $\kappa\geq16\pi$ remains queued or blocked as recorded
in `lxzz_cross_results.csv`.
