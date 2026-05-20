---
name: fem-verify
description: Create numerical experiments to verify FEM/DDM code — manufactured solutions, convergence rate tests, and benchmark comparisons.
allowed-tools: [Read, Write, Edit, Bash, Grep, Glob]
---

# FEM Verification Skill — Numerical Experiments

You write numerical experiments that verify the correctness and convergence of FEM/DDM implementations. Every piece of assembly code must have a corresponding verification.

## Verification Checklist

For every new FEM/DDM component, verify:
1. **Sparse pattern** — Is the matrix sparsity pattern correct?
2. **Patch test** — Does it pass the constant-strain patch test?
3. **Convergence rate** — Does error decay at the expected rate under uniform refinement?
4. **Manufactured solution** — For a known `u_exact`, does the discrete solution converge to it?

## Standard Test Template

```matlab
function verify_<ComponentName>()
    % VERIFY_<COMPONENTNAME>  Numerical verification of <component>.
    %
    % Tests:
    %   1. Sparsity pattern correctness
    %   2. Constant patch test
    %   3. Convergence rate (expected: O(h^k) in L2, O(h^{k-1}) in H1)

    %% ---- Test 1: Sparsity Pattern ----
    fprintf('Test 1: Sparsity pattern... ');
    mesh = squaremesh([0 1 0 1], 0.25);
    A = assembleStiffness(mesh, pde_poisson());
    % Verify: symmetric, positive semi-definite, correct bandwidth
    assert(issymmetric(A), 'Stiffness matrix must be symmetric');
    fprintf('PASSED\n');

    %% ---- Test 2: Patch Test ----
    fprintf('Test 2: Patch test... ');
    % ... verify that a linear solution is reproduced exactly
    fprintf('PASSED\n');

    %% ---- Test 3: Convergence Rate ----
    fprintf('Test 3: Convergence rate...\n');
    h = zeros(1,4);  errL2 = zeros(1,4);  errH1 = zeros(1,4);
    for k = 1:4
        mesh = squaremesh([0 1 0 1], 2^(-k-1));
        uh = femSolver(mesh, pde_manufactured());
        [errL2(k), errH1(k)] = computeError(mesh, uh, @u_exact, @grad_u_exact);
        h(k) = max(mesh.area)^(1/2);
        if k > 1
            rateL2 = log(errL2(k)/errL2(k-1)) / log(h(k)/h(k-1));
            fprintf('  h=%.4e  L2-err=%.4e  rate=%.2f\n', h(k), errL2(k), rateL2);
        end
    end
    assert(rateL2 > 1.8, 'Expected O(h^2) convergence in L2, got %.2f', rateL2);
    fprintf('Test 3: PASSED\n');
end
```

## Reusable Test Utilities

Keep these in a shared `tests/` or `verify/` directory:

| Utility | Purpose |
|----------|---------|
| `squaremesh(bbox, h)` | Uniform square mesh (parameterized by mesh size) |
| `pde_poisson()` | Standard Poisson PDE struct `{-Delta u = 1, u|_d = 0}` |
| `pde_manufactured()` | PDE with known exact solution (e.g., `u = sin(pi*x)*sin(pi*y)`) |
| `computeError(mesh, uh, u_exact, grad_exact)` | Compute L2 and H1 error norms via quadrature |
| `computeConvergenceRate(h, err)` | Fit convergence rate from mesh-size vs error data |

## Convergence Rate Reference

| Element | Degree | L2 rate | H1 rate |
|---------|--------|---------|---------|
| P1 (linear) | k=1 | O(h²) | O(h) |
| P2 (quadratic) | k=2 | O(h³) | O(h²) |
| P3 (cubic) | k=3 | O(h⁴) | O(h³) |

## For DDM Verification

- **Two-subdomain test:** Known interface solution, verify Schur complement
- **Multi-domain patch test:** Constant solution across subdomains
- **Iteration count:** Number of Schwarz iterations should be bounded independent of mesh size

## Helmholtz ORAS Verification Rules

For ORAS Helmholtz experiments, use linear nodal partition-of-unity weights by default. In this codebase, attach them with:

```matlab
parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);
parts = linearPartitionOfUnity2D(parts, bbox, gridSize, delta);
```

Use equal `1/nodeCount` weights only for an explicit comparison or when the user asks for plain RAS weighting. The Gander-Gong-Graham-Spence variational ORAS operator relies on weighted prolongation by a nodal partition of unity.

The POU must be genuinely piecewise linear across the full overlap. For neighbouring strips with overlap extension `delta` on each side, the two weights on `[x_interface-delta, x_interface+delta]` should be linear and sum to one. Do not use plateau weights that are later normalized; that produces a non-linear rational transition and does not match the GGS setup.

For Helmholtz iteration tables, explicitly enforce the resolution condition. When reproducing GGS Tables 5-8, sweep

```matlab
h = 2*pi/(q*k),  q in [10, 20, 40, 80]
```

or document clearly when a smaller compact run skips large `q` because MATLAB dense `E` or high-order sparse solves become too expensive. Do not replace this sweep by a single "about 10 points per wavelength" mesh unless the user requests a quick smoke test.

For power-norm studies, record both `||E||` and `||E^N||` in the `k`-weighted H1 norm induced by `K + k^2*M`, and report the partition shape, overlap convention, polynomial degree, `k`, and `h`.

For strip GGS iteration experiments, interpret overlap width as the total overlap between adjacent original subdomains. Thus "overlap width 1/2" means extension `1/4` on each side of each non-overlapped strip.

For GGS reproduction experiments, make the artificial subdomain boundaries and extended boundaries mesh-resolved. The paper assumes each `partial Omega_j` is resolved by the finite-element mesh. For the 8-strip domain `(0,16/3) x (0,1)` with `H=2/3` and extension `1/4`, use a mesh size of the form `h = 1/(12*m)` near the target `2*pi/(q*k)`, so both `H/h` and `(1/4)/h` are integers. For checkerboards with `gridN` subdomains per side and overlap `H/4`, use `h = 1/(4*gridN*m)` near the target.

Be explicit about the strip overlap convention. The user convention "overlap width 1/2" means extension `1/4` on each side; in diagnostics on this codebase, extension `1/2` on each side gives iteration counts much closer to Table 5. Use `ORAS_LARGEK_STRIP_OVERLAP_EXTENSION` to test both and report which convention is used.

For large wave-number Helmholtz runs, especially `k > 100`, estimate memory before running. Include global DOFs, polynomial degree, partition shape, approximate sparse matrix storage, local LU storage, and GMRES Krylov storage. Do not form dense `E` at large `k`; restrict large cases to Richardson/GMRES iteration tables unless the dense `E` memory estimate is clearly safe.

For large-k ORAS iteration batches, use `verify/verify_oras_largek_iterations.m`. It checkpoints after every completed case and supports environment switches:

```matlab
setenv('ORAS_LARGEK_KVALS','40 80 120');
setenv('ORAS_LARGEK_DEGREES','1 2 3');
setenv('ORAS_LARGEK_QVALS','10');
setenv('ORAS_LARGEK_PARPOOL','off');     % off, auto, or on
setenv('ORAS_LARGEK_WORKERS','4');       % optional
setenv('ORAS_LARGEK_TAG','run1');        % avoids checkpoint overwrites
```

Parallel MATLAB batches are allowed for these long experiments when CPU is available. Use distinct `ORAS_LARGEK_TAG` values so concurrent runs write separate `.mat` and `.md` checkpoint files. Keep a watchdog or command timeout for every long batch; for the current large-k studies, use a 6-hour wall-time limit unless the user gives a different limit.

When partial checkpoint results already make the numerical conclusion clear, stop and report the partial table to the user instead of waiting for the full sweep. This is especially important if the results closely match the target Gander-Gong-Graham-Spence tables, or if they show a systematic distinction that needs user advice before spending more CPU time.
