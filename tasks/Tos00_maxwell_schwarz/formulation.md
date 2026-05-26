Reproduction target: Tables 1-4 in Toselli, *Overlapping Schwarz methods for Maxwell's equations in three dimensions*.

Created: 2026-05-26
Updated: 2026-05-26
Verification entry point: `verify/verify_toselli_maxwell_schwarz.m`
Main utilities: `assembleMaxwell3D`, `nedelecAdditiveSchwarz3D`, `nedelecTwoLevelASM3D`, `prolongateNestedNed1`, `pcgLanczosCondition`

# Positive Definite Maxwell Schwarz Reproduction

## PDE Statement

The reproduced discrete problem is Toselli's Dirichlet positive definite Maxwell problem on the unit cube:

Find `u in V_h subset H_0(curl, Omega)` such that

```text
a(u,v) = (f,v),  for all v in V_h,
a(u,v) = eta1 * (u,v)_L2 + eta2 * (curl u, curl v)_L2.
```

The boundary condition is the perfect conductor condition `u x n = 0` on `partial Omega`, imposed by removing boundary edge DOFs.

## Discrete System

Use lowest-order first-kind Nedelec tetrahedral elements already implemented in the repo. The matrix is

```text
A = eta1 * M + eta2 * C,
M_ij = int_Omega phi_i . phi_j dx,
C_ij = int_Omega curl(phi_i) . curl(phi_j) dx.
```

After boundary edge elimination, PCG solves `A_ff u_f = b_f`.

## Schwarz Operators

The one-level additive Schwarz preconditioner is

```text
M_AS^{-1} = sum_i R_i^T A_i^{-1} R_i.
```

For Nedelec edge elements, the local interior DOF set contains the global free edges whose full tetrahedral support is inside the overlapping subdomain.

The two-level preconditioner is

```text
M_2AS^{-1} = P_0 A_0^{-1} P_0^T + sum_i R_i^T A_i^{-1} R_i,
A_0 = P_0^T A_ff P_0.
```

The coarse space uses NE1 on the `m^3` coarse cube mesh. `P_0` is the nested NE1 interpolation matrix, computed by exact edge moments of coarse basis functions over each fine edge.

## Table Parameters

- Tables 1-2: `eta1=1`, `eta2=1`; rows are `(n,m)` with fine mesh `n^3` cubes and `m^3` cubical subdomains; columns are `H/delta in {8,4,2,4/3}`.
- Table 1: one-level additive Schwarz.
- Table 2: two-level additive Schwarz.
- Tables 3-4: fixed `n=16`, `eta2=1`; vary `eta1`; compare `m=2` and `m=4` rows.
- A table cell is valid only when `delta/h = n/(m*(H/delta))` is an integer.

## Current Implementation Notes

The reproduction driver uses a deterministic synthetic right-hand side because the article reports CG iteration counts but does not specify the right-hand side. Condition numbers estimated from CG/Lanczos data are the primary comparison target.
