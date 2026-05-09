---
name: fem-assemble
description: Write reusable, vectorized FEM assembly routines (stiffness, mass, load) following iFEM patterns. Ensures consistent API across all element types.
allowed-tools: [Read, Write, Edit, Bash, Grep, Glob, WebFetch, WebSearch]
---

# FEM Assembly Skill — Vectorized & Reusable

You are an expert in writing vectorized MATLAB FEM assembly code following the patterns established in Chen Long's iFEM package. Your goal is to produce code that is reusable across element types, problem dimensions, and quadrature rules.

## Code Reuse Principles

1. **One function, one responsibility** — Stiffness, mass, load, and BC application each get their own function. Never combine them.
2. **Uniform API** — Every assembly function takes `(mesh, pde)` structs and returns a sparse matrix or vector:
   ```matlab
   function A = assembleStiffness(mesh, pde)
   function M = assembleMass(mesh, pde)
   function b = assembleLoad(mesh, pde)
   function [A, b] = applyDirichletBC(A, b, mesh, pde)
   ```
3. **mesh struct** — All geometric data in one struct: `node`, `elem`, `bdFlag`, `area`, etc.
4. **pde struct** — All PDE data in one struct: `coefficient functions`, `boundary conditions`, `source term`.
5. **No hardcoded constants** — Quadrature order, element type, and material parameters come from `mesh` or `pde`.

## Vectorized Assembly Pattern (2D Triangular Elements)

```matlab
function A = assembleStiffness(mesh, pde)
    % Input structs
    % mesh.node: N x 2 vertex coordinates
    % mesh.elem: NT x 3 element connectivity (1-indexed)
    % pde.coef:  coefficient function handle @(x,y) or scalar

    N = size(mesh.node, 1);
    NT = size(mesh.elem, 1);

    % ---- Quadrature (reusable — call a shared quadrature function) ----
    [lambda, weight] = quadpts(2);  % order-2 quadrature on triangle

    % ---- Geometry (vectorized across elements) ----
    v1 = mesh.node(mesh.elem(:,1), :);
    v2 = mesh.node(mesh.elem(:,2), :);
    v3 = mesh.node(mesh.elem(:,3), :);

    % ---- Preallocate for sparse assembly ----
    ii = zeros(9*NT, 1);
    jj = zeros(9*NT, 1);
    ss = zeros(9*NT, 1);

    % ---- Element-wise computation ----
    for p = 1:length(weight)
        % Evaluate basis gradients at quadrature point (vectorized)
        [Dphi, area] = gradbasis(v1, v2, v3);  % [#elem x 3 x 2]

        % Coefficient at quadrature point
        pxy = lambda(p,1)*v1 + lambda(p,2)*v2 + lambda(p,3)*v3;
        coef_val = pde.coef(pxy(:,1), pxy(:,2));  % [#elem x 1]

        % Local stiffness contribution: weight * coef * (Dphi * Dphi') * area
        % Vectorized over elements using Einstein summation pattern
        for i = 1:3
            for j = 1:3
                loc = Dphi(:,:,i) .* Dphi(:,:,j);  % [#elem x 2]
                s = weight(p) * coef_val .* sum(loc, 2) .* area;
                % store into ii, jj, ss at the right block
            end
        end
    end

    % ---- Global assembly (reusable pattern) ----
    A = sparse(ii, jj, ss, N, N);
end
```

## Key Reusable Utilities to Reference

Always use from the shared utility pool rather than re-implementing:

| Utility | Purpose |
|----------|---------|
| `quadpts(order)` | Gauss quadrature on reference element |
| `gradbasis(v1,v2,v3)` | Gradients of barycentric coordinates |
| `setboundary(mesh)` | Extract boundary edges |
| `sparse(ii,jj,ss,N,N)` | Global assembly (built-in, keep this pattern) |
| `applyDirichletBC(A,b,mesh,pde)` | Enforce Dirichlet conditions |

## When Writing Assembly Code

- **NEVER** loop over elements in MATLAB — vectorize across elements
- Use `repmat`, `reshape`, `bsxfun` (or implicit broadcasting) over `for` loops
- Pre-allocate the index arrays `ii, jj, ss` and fill in blocks
- Use `sparse` with accumulated duplicates for the final matrix
- Extract repeated geometric computations to shared utility functions
- Comment the mathematical formula being discretized in a one-liner above each block
