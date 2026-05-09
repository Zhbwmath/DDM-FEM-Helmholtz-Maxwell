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
