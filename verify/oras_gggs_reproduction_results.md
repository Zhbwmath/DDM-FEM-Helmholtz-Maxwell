# ORAS Gander-Gong-Graham-Spence Reproduction Notes

Date: 2026-05-19

References:
- Gong, Gander, Graham, Spence, "A variational interpretation of Restricted Additive Schwarz with impedance transmission condition for the Helmholtz problem", Table 1.
- Gong, Graham, Spence, "Convergence of ORAS for discrete Helmholtz problems", Tables 5-8.

## Implementation Rules Now Used

ORAS reproduction runs use linear nodal partition-of-unity weights by default via:

```matlab
parts = linearPartitionOfUnity2D(parts, bbox, gridSize, overlap);
```

Equal-count weights `1/nodeCount` are kept only as a fallback in `orasHelmholtz` when no `parts(s).weightFun` is present. They are not used for GGS reproduction because the variational ORAS operator uses weighted prolongation by a nodal partition of unity.

For iteration tables, the script now enforces the Helmholtz resolution sweep

```matlab
h = 2*pi/(q*k),  q in [10, 20, 40, 80].
```

Rows whose estimated DOF exceeds the compact MATLAB threshold are printed as `skip`; they are no longer silently replaced by a coarser "about 10 ppw" mesh.

## Power-Norm Reproduction

Norm: k-weighted H1 norm induced by `D_k = K + k^2 M`.

Strip domain `(0, 2N/3) x (0, 1)`, P2, overlap extension `1/6`:

| Nsub | k | 1/h | DOF | ||E|| | ||E^N|| | GMRES |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 10 | 14 | 1131 | 3.36 | 0.278 | 5 |
| 2 | 20 | 14 | 1131 | 3.56 | 0.296 | 6 |
| 4 | 10 | 14 | 2175 | 3.35 | 0.440 | 11 |
| 4 | 20 | 14 | 2175 | 3.81 | 0.492 | 14 |

Checkerboard unit square, P2, overlap `H/4`:

| grid | k | 1/h | DOF | ||E|| | ||E^N|| | GMRES |
|---:|---:|---:|---:|---:|---:|---:|
| 2x2 | 10 | 14 | 841 | 3.55 | 0.030 | 8 |
| 2x2 | 20 | 14 | 841 | 3.60 | 0.070 | 9 |
| 3x3 | 10 | 14 | 841 | 3.33 | 0.991 | 16 |
| 3x3 | 20 | 14 | 841 | 3.39 | 0.033 | 17 |

## Iteration Reproduction With Resolution Sweep

Residual reduction target: `1e-6` in Euclidean norm.

Strip domain `(0,16/3) x (0,1)`, 8 strips:

| k | p | q | 1/h | DOF | Richardson | GMRES |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1 | 10 | 16 | 1462 | 100 | 23 |
| 10 | 1 | 20 | 32 | 5676 | >120 | 23 |
| 10 | 2 | 10 | 16 | 5643 | 59 | 22 |
| 20 | 1 | 10 | 32 | 5676 | 111 | 39 |

Higher `q` rows and high-order strip rows exceed the compact MATLAB threshold. Example skipped estimates:

| k | p | q | estimated DOF |
|---:|---:|---:|---:|
| 10 | 1 | 40 | 22295 |
| 10 | 2 | 20 | 22295 |
| 20 | 1 | 20 | 22295 |
| 20 | 2 | 10 | 22295 |
| 20 | 3 | 80 | 3126046 |

Checkerboard unit square, `H ~ k^{-0.4}`, overlap `H/4`:

| k | p | q | grid | 1/h | DOF | Richardson | GMRES |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1 | 10 | 3x3 | 16 | 289 | 38 | 15 |
| 10 | 1 | 20 | 3x3 | 32 | 1089 | 38 | 16 |
| 10 | 1 | 40 | 3x3 | 64 | 4225 | 46 | 15 |
| 10 | 2 | 10 | 3x3 | 16 | 1089 | 36 | 16 |
| 10 | 2 | 20 | 3x3 | 32 | 4225 | 37 | 16 |
| 10 | 3 | 10 | 3x3 | 16 | 2401 | 78 | 17 |
| 10 | 3 | 20 | 3x3 | 32 | 9409 | >120 | 19 |
| 20 | 1 | 10 | 3x3 | 32 | 1089 | 16 | 14 |
| 20 | 1 | 20 | 3x3 | 64 | 4225 | 18 | 14 |
| 20 | 2 | 10 | 3x3 | 32 | 4225 | 18 | 16 |
| 20 | 3 | 10 | 3x3 | 32 | 9409 | 26 | 16 |

## 2026-05-20 Diagnosis Against GGS Discrete ORAS

The discrete ORAS paper defines

```text
B_h^{-1} = sum_j \tilde R_{h,j}^T A_{h,j}^{-1} R_{h,j},
u_h^{n+1} = u_h^n + B_h^{-1}(F_h - A_h u_h^n),
```

where `A_{h,j}` is the local Helmholtz impedance operator on a mesh-resolved subdomain and
`\tilde R_{h,j}^T v_j = R_{h,j}^T Pi_h(chi_j v_j)`.

Three code issues were identified and fixed:

1. `linearPartitionOfUnity2D` used plateau ramps normalized afterwards. The normalized weights were not genuinely linear across the full overlap. It now uses 1D linear ramps that sum to one across the whole overlap.
2. The large-k strip/checkerboard runner did not force artificial subdomain boundaries and extended boundaries to align with the structured mesh. The runner now chooses aligned mesh sizes for GGS-style cases.
3. `assembleBoundaryMass2D` had P2/P3 edge DOFs ordered inconsistently with the 1D Lagrange basis used for boundary quadrature. This corrupted high-order impedance terms on global and local ORAS boundaries.

Focused diagnostics after these fixes:

| convention | k | p | strip extension | 1/h | Richardson | GMRES |
|---|---:|---:|---:|---:|---:|---:|
| user total-overlap 1/2 | 20 | 1 | 1/4 | 36 | 25 | 24 |
| user total-overlap 1/2 | 40 | 1 | 1/4 | 60 | 25 | 24 |
| user total-overlap 1/2 | 80 | 1 | 1/4 | 132 | 26 | 24 |
| paper-count convention | 20 | 1 | 1/2 | 30 | 16 | 15 |
| paper-count convention | 40 | 1 | 1/2 | 66 | 17 | 16 |
| user total-overlap 1/2 | 20 | 2 | 1/4 | 36 | 22 | 20 |

The strip discrepancy is now mainly an overlap-convention issue. Extension `1/2` on each side gives iteration counts close to Table 5 (`14 (14)` at `k=20`, `14 (14)` at `k=40`), while extension `1/4` on each side gives stable but slower counts around `25 (24)`.

After vectorizing high-order mass/boundary-mass assembly and switching ORAS local solves to vector-form sparse LU with row/column permutations, the P2/P3 strip cases also match the fast ORAS/GMRES behavior:

| k | p | q | shape | 1/h | DOF | Richardson | GMRES |
|---:|---:|---:|---|---:|---:|---:|---:|
| 20 | 2 | 10 | strip | 30 | 19581 | 16 | 15 |
| 20 | 2 | 10 | grid 3x3 | 36 | 5329 | 18 | 14 |
| 20 | 3 | 10 | strip | 30 | 43771 | 16 | 15 |
| 20 | 3 | 10 | grid 3x3 | 36 | 11881 | 17 | 14 |
| 40 | 2 | 10 | strip | 66 | 93765 | 15 | 16 |
| 40 | 2 | 10 | grid 4x4 | 64 | 16641 | 19 | 19 |
| 40 | 3 | 10 | strip | 66 | 210343 | 15 | 16 |
| 40 | 3 | 10 | grid 4x4 | 64 | 37249 | 19 | 19 |

## Current Distinctions From Published Tables 5-8

1. The previous script did not genuinely use the Table 5-style resolution sweep. It now does.
2. With linear POU, the power-norm behavior qualitatively matches the papers: `||E|| > 1` but `||E^N|| < 1`.
3. Richardson iteration counts are still larger than the published tables in the strip cases.
4. Many published-resolution strip cases are too large for this compact MATLAB harness. Running them faithfully will require accepting much longer sparse-only runs, avoiding explicit dense `E`, or moving the experiment to a more scalable backend.
5. The strip overlap convention still needs confirmation: the script currently treats overlap width `1/2` as total adjacent overlap, i.e. extension `1/4` on each side.
