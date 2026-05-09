---
name: fem-component
description: Create new reusable FEM/DDM components (element types, quadrature rules, solvers, preconditioners) following the project's established patterns and API conventions.
allowed-tools: [Read, Write, Edit, Bash, Grep, Glob, WebFetch, WebSearch]
---

# FEM Component Skill — Reusable Building Blocks

You create new FEM/DDM components that integrate with the existing codebase. Every component follows the project's struct-based API conventions and is designed for composition.

## Component Categories

### 1. Mesh Components
```matlab
function mesh = squaremesh(bbox, h)
function mesh = circlemesh(center, radius, h)
function mesh = uniformrefine(mesh)        % h-refinement
function mesh = bisect(mesh, markedElem)   % adaptive refinement
```

**mesh struct fields:** `node` (Nx2), `elem` (NTx3 or NTx4), `bdFlag` (NTx3), `area` (NTx1), `edge` (NE x2), `neighbor` (NTx3)

### 2. Quadrature Components
```matlab
function [lambda, weight] = quadpts(order)       % Triangle reference
function [lambda, weight] = quadpts1d(order)     % 1D reference [-1,1]
function [xi, weight] = quadptsCube(order, dim)  % Hypercube reference
```

**Convention:** Outputs always `(points, weights)` where points is Nq x (dim+1) for barycentric, Nq x dim for Cartesian.

### 3. Finite Element Space Components
```matlab
function [phi, Dphi, volume] = basisP1(node, elem, quadpts)
function [phi, Dphi, volume] = basisP2(node, elem, quadpts)
function [phi, Dphi, volume] = basisCR(node, elem, quadpts)   % Crouzeix-Raviart
```

**Return convention:**
- `phi`:  `NT x Nq x Nlb` — basis values at quadrature points
- `Dphi`: `NT x Nq x Nlb x dim` — basis gradients
- `volume`: `NT x 1` — element volumes/areas

### 4. PDE Components (pde struct)
```matlab
function pde = pde_poisson()
function pde = pde_elasticity(mu, lambda)  % Lame parameters
function pde = pde_advection_diffusion(epsilon, beta)
```

**pde struct fields:** `type` (string), `coef` (@(x,y) or scalar), `source` (@(x,y)), `dirichlet` (struct with `node` and `value`), `neumann` (struct with `edge` and `value`), `exact` (@(x,y) optional, for verification).

### 5. Solver Components
```matlab
function uh = femSolver(mesh, pde)         % Direct solve (A\b)
function uh = femSolverCG(mesh, pde, tol)  % Conjugate gradient
function uh = ddmSchur(mesh, pde, subdomains)  % DDM with Schur complement
```

### 6. DDM Components
```matlab
function [submesh, interface] = partitionMesh(mesh, nSubdomains)
function S = assembleSchurComplement(submesh, interface, pde)
function uh = ddmSchwarz(submesh, interface, pde, maxIter, tol)
```

## Code Reuse Rules When Creating Components

1. **Check for existing utilities first** — Use `Grep` to search for `quadpts`, `gradbasis`, `assembleStiffness` etc. before writing new ones.
2. **Match the existing API** — If `assembleStiffness` already takes `(mesh, pde)`, your new component must too.
3. **Put shared math in utilities, not in assembly** — `gradbasis` should be a standalone function, not embedded inside `assembleStiffness`.
4. **Struct over positional arguments** — If a function needs 5+ arguments, wrap them in a struct.
5. **Comment the discretization formula** — One line above the computation showing the mathematical formula.
6. **Write the verification in the same commit** — Every component ships with its `verify_*.m`.

## Before Writing Any New Component

Search the existing codebase for:
- Similar element types (can they be generalized instead?)
- Existing quadrature rules (is the order already implemented?)
- Existing PDE structs (can `pde.coef` be replaced with a function handle rather than creating a new struct type?)

The goal is **fewer files, more reusable functions.**
