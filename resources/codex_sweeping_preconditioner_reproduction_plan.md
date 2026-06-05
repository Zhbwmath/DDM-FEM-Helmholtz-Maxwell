# Codex Reproduction Plan: Sweeping Preconditioners for the Helmholtz Equation

## Objective

Reproduce the main numerical experiments from the following two papers:

1. B. Engquist and L. Ying, *Sweeping Preconditioner for the Helmholtz Equation: Moving Perfectly Matched Layers*.
2. B. Engquist and L. Ying, *Sweeping Preconditioner for the Helmholtz Equation: Hierarchical Matrix Representation*.

The priority is to implement a correct and testable 2D moving-PML sweeping preconditioner first. Then add a simplified hierarchical-matrix version, preferably HODLR/H-matrix rank-compressed Schur complements. Exact wall-clock times from the papers are not required, because hardware, programming language, and sparse/direct-solver libraries differ. The main reproducibility targets are:

- correct Helmholtz/PML discretization;
- correct layer-wise sweeping preconditioner structure;
- GMRES iteration counts that are roughly independent of frequency and problem size;
- setup and solve time scaling close to the paper trends;
- tables analogous to the paper tables.

---

## Non-negotiable implementation principles

Do not reproduce the experiment by simply calling a global sparse direct solve as the solver. A global direct solve may only be used for reference validation on small grids.

Implement the preconditioner as a `LinearOperator` or equivalent object whose `apply(rhs)` executes an approximate block inverse/sweeping solve.

The implementation must separate:

1. global Helmholtz matrix assembly;
2. local moving-PML strip/slab construction;
3. preconditioner setup;
4. preconditioner application;
5. GMRES experiment driver;
6. table/plot generation.

Use deterministic random seeds for random velocity fields and randomized low-rank compression.

---

## Recommended repository structure

```text
sweeping-preconditioner-repro/
  README.md
  requirements.txt
  src/
    grid.py
    pml.py
    velocities.py
    sources.py
    assemble_helmholtz.py
    layer_blocks.py
    moving_pml_preconditioner.py
    hmatrix_preconditioner.py
    gmres_driver.py
    diagnostics.py
    plotting.py
  experiments/
    run_2d_moving_pml_vary_omega.py
    run_2d_moving_pml_vary_q.py
    run_2d_hmatrix_pilot.py
    run_3d_moving_pml_pilot.py
  results/
    tables/
    figures/
  tests/
    test_assembly_small.py
    test_layer_blocks.py
    test_moving_pml_action.py
    test_gmres_reference.py
```

Suggested Python stack:

```text
numpy
scipy
matplotlib
pandas
pytest
```

If using MATLAB instead, keep the same module decomposition in separate `.m` files.

---

## Mathematical problem

Solve the variable-coefficient Helmholtz equation on the unit square or unit cube:

```math
\Delta u(x) + \frac{\omega^2}{c(x)^2}u(x)=f(x).
```

Use exterior PML/absorbing layers on the boundary. Use a 5-point finite-difference stencil in 2D and a 7-point finite-difference stencil in 3D.

The paper convention is:

```math
\omega/(2\pi) = \text{number of wavelengths across the unit domain},
```

and the number of grid points per wavelength is `q`. Therefore use

```text
n = q * (omega / (2*pi)) - 1
h = 1 / (n + 1)
N = n^2 in 2D, N = n^3 in 3D.
```

For example, in 2D with `q=8` and `omega/(2*pi)=16`, use `n=127`.

---

## PML model to implement

Use complex coordinate stretching. Define

```math
s_j(t)=\left(1+i\frac{\sigma_j(t)}{\omega}\right)^{-1}.
```

Use a quadratic absorption profile near the boundary:

```math
\sigma(t)=
\begin{cases}
(C/\eta)((t-\eta)/\eta)^2, & 0\le t\le \eta,\\
0, & \eta<t<1-\eta,\\
(C/\eta)((t-1+\eta)/\eta)^2, & 1-\eta\le t\le 1.
\end{cases}
```

For one-sided local moving PML strips, use the same quadratic profile only on the artificial absorbing side.

Start with:

```text
C = 4 or 6
eta = pml_width_grid * h
pml_width_grid = b
```

For the 2D moving-PML reproduction, use:

```text
alpha = 2
b = 12
q = 8
gmres relative residual tolerance = 1e-3
```

For the 3D moving-PML pilot, use:

```text
alpha = 1
b = 6
d = 3
gmres relative residual tolerance = 1e-3
```

The preconditioner should be constructed for the damped operator

```math
\Delta u + \frac{(\omega+i\alpha)^2}{c(x)^2}u=f,
```

but GMRES should solve the original undamped system.

---

## Velocity fields for the 2D experiments

Implement three velocity fields on `(0,1)^2`.

### Velocity field 1: Gaussian converging lens

Use center `(r1,r2)=(1/2,1/2)`:

```math
c(x_1,x_2)=\frac{4}{3}\left(1-\frac12\exp(-32((x_1-r_1)^2+(x_2-r_2)^2))\right).
```

### Velocity field 2: vertical Gaussian waveguide

```math
c(x_1,x_2)=\frac{4}{3}\left(1-\frac12\exp(-32(x_1-1/2)^2)\right).
```

### Velocity field 3: random smooth velocity field

Generate values in `(0.7,1.3)` with correlation length approximately `1/16` in 2D. Practical implementation:

1. draw iid standard normal noise on the grid;
2. smooth with a Gaussian filter of physical width approximately `1/16`;
3. rescale linearly to `[0.7,1.3]`;
4. fix the seed, for example `seed=202601`.

This will not exactly match the paper’s random field, but it should reproduce the qualitative test.

---

## Source terms for the 2D experiments

### Source 1: narrow Gaussian point source

Center `(r1,r2)=(1/2,1/8)`:

```math
f_1(x_1,x_2)=\exp\left(-(4\omega/\pi)^2((x_1-r_1)^2+(x_2-r_2)^2)\right).
```

### Source 2: Gaussian wave packet

Center `(r1,r2)=(1/8,1/8)` and direction

```math
(d_1,d_2)=\frac{1}{\sqrt{2}}(1,1).
```

Use

```math
f_2(x_1,x_2)=\exp(-4\omega((x_1-r_1)^2+(x_2-r_2)^2))
\exp(i\omega(x_1d_1+x_2d_2)).
```

Normalize the source vector if needed, but report clearly whether normalization is used.

---

## Global matrix assembly in 2D

Implement a centered finite-difference PML Helmholtz matrix. Use row-major layer ordering:

```text
layer m = all grid points (i,m), i=1,...,n
unknown index = i + n*m, with zero-based m
```

Construct sparse matrix `A` for the original operator with frequency `omega`.
Construct sparse matrix `A_alpha` for the preconditioner with frequency `omega + 1j*alpha`.

For debugging, also implement a no-PML version of the matrix and verify the stencil on very small grids.

Acceptance tests:

- matrix shape is `(n*n, n*n)`;
- each interior row has at most 5 nonzeros;
- PML matrix is complex symmetric up to numerical tolerance if using the symmetric PML form;
- block extraction produces block tridiagonal structure in the sweep direction.

---

## Layer block extraction

Implement a function:

```python
def extract_layer_blocks(A, n):
    """Return diagonal blocks A_mm and off-diagonal blocks A_m_mplus1, A_mplus1_m."""
```

For the 5-point stencil and row-wise layers:

- `A_mm` is an `n x n` tridiagonal sparse matrix;
- `A_m,m+1` and `A_m+1,m` are diagonal sparse matrices;
- non-neighbor layer blocks should be zero.

Add a test that reconstructs `A` from layer blocks and checks the reconstruction error.

---

## Moving-PML sweeping preconditioner: 2D primary implementation

Implement the preconditioner as an approximate block inverse of `A_alpha`.

### Local operator \tilde{T}_m

For each layer `m`, define a local strip containing the current layer and `b` layers behind it:

```text
rows = max(0, m-b+1), ..., m
```

Build a local Helmholtz/PML matrix `H_m` on this strip with an artificial PML on the lower/back side of the strip. Factorize `H_m` once during setup.

The approximate Schur-complement inverse action is:

```python
def apply_T_tilde(m, g_m):
    # g_m has length n and is supported on the current layer m
    local_rhs = zeros(local_strip_size)
    insert g_m on the last/current strip layer
    local_sol = solve_factorized_H_m(local_rhs)
    return extract local_sol on the current strip layer
```

Use sparse LU (`scipy.sparse.linalg.splu`) for the first implementation. In 2D the strip is quasi-1D, so this is sufficient for correctness. Later, use banded LU if desired.

### Preconditioner setup

```python
def setup_moving_pml_preconditioner(A_alpha, omega, alpha, c, n, b):
    for m in range(n):
        build H_m
        factorize H_m
        store factorization and extraction maps
```

Record setup time.

### Preconditioner application

Start with the simpler one-layer version (`d=1`):

```python
def apply_M(rhs):
    u = split rhs into n layers

    # Forward sweep / approximate elimination
    for m in range(0, n-1):
        tmp = apply_T_tilde(m, u[m])
        u[m+1] -= A[m+1,m] @ tmp

    # Diagonal layer solves
    for m in range(0, n):
        u[m] = apply_T_tilde(m, u[m])

    # Backward substitution
    for m in range(n-2, -1, -1):
        tmp_rhs = A[m,m+1] @ u[m+1]
        u[m] -= apply_T_tilde(m, tmp_rhs)

    return concatenate u
```

After this works, implement grouped layers `d=12` for closer reproduction of the paper.

---

## GMRES driver

Use left preconditioning:

```math
M_\alpha A u = M_\alpha f.
```

In SciPy this can be implemented by giving `M=M_alpha` to `scipy.sparse.linalg.gmres`.

Record:

```text
omega_over_2pi
q
n
N
velocity_id
source_id
Tsetup
Niter
Tsolve
relative_residual
```

Use relative residual tolerance `1e-3`. Also compute and store the true residual

```math
\|Au-f\|_2/\|f\|_2
```

after GMRES terminates.

---

## Main 2D moving-PML experiments to reproduce

### Experiment A: frequency scaling, analogous to Tables 1--3 in the moving-PML paper

Use:

```text
velocity_id in {1,2,3}
source_id in {1,2}
omega_over_2pi in {16, 32, 64, 128, 256}
q = 8
alpha = 2
b = 12
d = 1 initially, d = 12 later
GMRES tolerance = 1e-3
```

For memory reasons, first run only up to `omega_over_2pi=64` or `128`. Then run `256` only after correctness and memory use are stable.

Expected qualitative result:

- moving-PML preconditioner: iteration counts roughly in the range 14--23 for the 2D tests;
- iteration counts should not grow significantly with `n`;
- setup time should scale roughly linearly with `N` in 2D, up to implementation overhead.

### Experiment B: points-per-wavelength scaling, analogous to Tables 4--6 in the moving-PML paper

Use:

```text
velocity_id in {1,2,3}
source_id in {1,2}
omega_over_2pi = 32
q in {8, 16, 32, 64}
alpha = 2
b = 12
GMRES tolerance = 1e-3
```

Expected qualitative result:

- iteration counts should be almost constant or grow slowly;
- solve time should scale near linearly with `N`, subject to the cost of the local sparse LU solves.

---

## 3D moving-PML pilot experiment

Do this only after 2D works.

Use:

```text
D = (0,1)^3
5? No: use 7-point finite-difference stencil
velocity_id in {1,2,3}
source_id in {1,2}
omega_over_2pi in {5,10,20}
q = 8
n in {39,79,159}
alpha = 1
b = 6
d = 3
GMRES tolerance = 1e-3
```

Sources:

- point source centered at `(1/2,1/2,1/4)`;
- Gaussian wave packet centered at `(1/2,1/4,1/4)` and pointing in direction `(0,1,1)/sqrt(2)`.

For local slab solves, use sparse LU first. A full multifrontal implementation is not required for the first reproduction. If `n=159` is too large, run only `n=39` and `n=79` and report the limitation.

Expected qualitative result:

- iteration counts around 11--15 for moving PML in 3D;
- setup much more expensive than 2D but still far below a global 3D direct solve.

---

## Hierarchical-matrix reproduction plan

The full hierarchical-matrix version is harder than the moving-PML version. Implement it in stages.

### H-matrix Stage 1: exact dense Schur complement on small grids

For small `n`, compute the exact block Schur recursion:

```math
S_1=A_{1,1},
\quad
T_1=S_1^{-1},
```

```math
S_m=A_{m,m}-A_{m,m-1}T_{m-1}A_{m-1,m},
\quad
T_m=S_m^{-1}.
```

Use this exact dense version only for small grids, for example `n <= 127`, to validate the block factorization.

### H-matrix Stage 2: HODLR compression of each \tilde{T}_m

Implement a simple one-dimensional binary cluster tree for the indices of one layer.

For a dense layer matrix `T_m`, store:

- diagonal blocks densely;
- off-diagonal blocks in truncated SVD form with fixed rank `r`.

Start with:

```text
r = 2 for 2D tests
leaf_size = 32 or 64
```

This will not be a full H-matrix implementation, but it is enough to test the main insight: off-diagonal layer Green's interactions are compressible.

### H-matrix Stage 3: approximate Schur recursion

Replace dense operations in the Schur recursion by compressed operations:

```math
\tilde{S}_m \approx A_{m,m}-A_{m,m-1}\tilde{T}_{m-1}A_{m-1,m},
\quad
\tilde{T}_m \approx \tilde{S}_m^{-1}.
```

If exact H-matrix inversion is too time-consuming to implement, use this fallback:

1. form `S_m` densely for moderate `n`;
2. invert it densely;
3. compress `T_m` into HODLR form;
4. use the compressed apply in the sweep.

This fallback will not reproduce linear-complexity setup, but it can reproduce the iteration-count behavior and the compression phenomenon.

### H-matrix experiments

Run the same 2D velocity/source tests as the moving-PML version:

```text
omega_over_2pi in {16,32,64,128,256}
q = 8
rank r = 2
GMRES tolerance = 1e-3
```

Expected qualitative result from the paper:

- H-matrix sweeping uses very small ranks in 2D;
- GMRES iteration counts are often around 2--5;
- setup can be more expensive than moving PML, but the preconditioner is usually more accurate.

If the full `omega_over_2pi=256` case is impossible, report results for the largest feasible case and include a scaling plot.

---

## Validation tests before large experiments

### Test 1: matrix assembly

For `n=8`, assemble `A` and verify:

- dimensions;
- sparsity pattern;
- stencil signs;
- block tridiagonal structure.

### Test 2: local strip solver

For `n=32`, `b=6`, choose a random layer vector `g_m`. Apply `T_tilde_m(g_m)`. Verify:

- output has length `n`;
- no NaNs/Infs;
- result changes smoothly when `b` increases.

### Test 3: exact Schur complement comparison

For `n=16` or `32`, compute exact `T_m` from the dense Schur recursion. Compare with moving-PML `T_tilde_m`:

```math
\frac{\|T_m g - \tilde{T}_m g\|_2}{\|T_m g\|_2}.
```

The error need not be tiny, but it should decrease when `b` increases or `alpha` is tuned.

### Test 4: preconditioned GMRES vs unpreconditioned GMRES

For `n=63`, compare:

- unpreconditioned GMRES iteration count;
- moving-PML preconditioned GMRES iteration count.

The preconditioned version should be dramatically better.

---

## Output format

Each experiment script should produce:

```text
results/tables/<experiment_name>.csv
results/tables/<experiment_name>.md
results/figures/<experiment_name>_velocity.png
results/figures/<experiment_name>_solution_source1.png
results/figures/<experiment_name>_solution_source2.png
results/figures/<experiment_name>_iterations_vs_n.png
results/figures/<experiment_name>_time_vs_N.png
```

Each Markdown table should contain:

```text
omega_over_2pi | q | n | N | velocity_id | source_id | Tsetup | Niter | Tsolve | final_true_residual
```

For comparison with the paper, also generate compact tables grouped by `velocity_id`, with columns:

```text
omega_over_2pi | q | N=n^2 | Tsetup | Niter source 1 | Tsolve source 1 | Niter source 2 | Tsolve source 2
```

---

## Performance and memory notes

First implementation can use sparse LU for every local strip. This is not optimal but is acceptable for reproducing the algorithmic structure.

Potential bottlenecks:

- storing one LU factorization per layer;
- Python overhead during GMRES if `apply_M` is not vectorized enough;
- repeated construction of local PML matrices.

Optimizations after correctness:

1. process multiple physical layers per sweep step (`d=12` in 2D moving-PML paper);
2. cache repeated local matrix sparsity patterns;
3. use banded LU for 2D local strips;
4. use compiled sparse solvers if available;
5. for 3D, use a multifrontal solver if available.

---

## Minimal milestone sequence

### Milestone 1: global Helmholtz/PML solver

- Assemble `A` and `A_alpha` in 2D.
- Solve one small problem using direct solve for reference.
- Plot velocity and solution for one source.

### Milestone 2: moving-PML local strip solvers

- Build and factorize `H_m` for all layers.
- Test `apply_T_tilde(m,g)`.

### Milestone 3: moving-PML preconditioner

- Implement `apply_M`.
- Use it inside GMRES.
- Run a small test: `omega_over_2pi=16`, `q=8`, velocity field 1, both sources.

### Milestone 4: reproduce 2D moving-PML tables

- Run velocity fields 1--3.
- Run frequency scaling and q-scaling experiments.
- Export tables and plots.

### Milestone 5: H-matrix pilot

- Implement dense Schur recursion for small `n`.
- Add HODLR compression of layer inverses.
- Demonstrate low-rank compression and reduced GMRES iterations.

### Milestone 6: 3D pilot

- Implement 7-point 3D matrix assembly.
- Implement moving-PML slab preconditioner with sparse LU.
- Run only feasible grid sizes.

---

## Expected numerical checkpoints

For 2D moving PML with `q=8`, `alpha=2`, `b=12`:

- frequency scaling uses `omega/(2*pi)=16,32,64,128,256` and `n=127,255,511,1023,2047`;
- expected GMRES iterations are roughly 14--23;
- iteration count should not grow proportionally to `n`.

For 2D H-matrix sweeping with rank around `r=2`:

- expected GMRES iterations are roughly 2--5 in the paper;
- exact timings are not a target unless an optimized H-matrix implementation is used.

For 3D moving PML:

- use `omega/(2*pi)=5,10,20`, `q=8`, `n=39,79,159`;
- expected GMRES iterations are roughly 11--15;
- if sparse LU memory prevents `n=159`, stop at `n=79` and report memory use.

---

## Reporting requirements

At the end, write `results/REPORT.md` with:

1. exact commands used;
2. software versions and hardware information;
3. all parameter values: `omega`, `q`, `n`, `b`, `alpha`, PML profile, GMRES tolerance;
4. tables of iteration counts and timings;
5. plots of iteration counts versus `n` and timings versus `N`;
6. a section explaining any deviations from the papers;
7. a section comparing moving-PML versus H-matrix preconditioners.

Use this interpretation for the final comparison:

- moving PML is easier to implement and more robust in 3D;
- H-matrix Schur complement approximation is more accurate in 2D but substantially harder to implement;
- both rely on the same layer-wise approximate block factorization and the half-space Green's function interpretation of Schur complement inverses.
