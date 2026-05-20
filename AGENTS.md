# DDM-FEM-Helmholtz-Maxwell Project

## MATLAB Execution

- **Always run MATLAB silently**: use `-nosplash -nodesktop -batch` — no console windows.
- **Never use `-noFigureWindows`** unless the user explicitly asks to suppress all graphics. Experiment figures must display.
- For running a script: `matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('script.m');"`
- For running inline code: `matlab -nosplash -nodesktop -batch "addpath(genpath('.')); <code>;"`
- Check for MATLAB at `C:\Program Files\MATLAB\R2024b\bin\matlab.exe` first, then fall back to `matlab` on PATH.

## HPC Rules (Workstation: 48 cores, 549 GB RAM)

- **Ask permission before any test estimated to use >100 GB memory.** Estimate memory usage and present it before running.
- **Use `parfor` for subdomain setup.** Start parpool before large runs: `parpool('local', feature('numcores'))`.
- **Vectorize all assembly.** Use iFEM-style one-shot `sparse(ii, jj, ss, N, N)`.
- **Avoid `containers.Map` with string keys for large meshes.** Use `ismember` on sorted edge lists or sparse hashing.
- **Memory estimation rule of thumb** for 2D Helmholtz with N nodes, NT elements, nSub subdomains:
  - Global sparse matrix A: ~7N nonzeros × 16 bytes (complex) ≈ 112N bytes
  - Edge list: ~6NT entries × 8 bytes ≈ 48NT bytes
  - Per subdomain (N/nSub nodes): ~3 × 7×N/nSub nonzeros × 8 bytes (real) + LU ≈ 500×(N/nSub)^1.5 bytes
  - Total ≈ 112N + 48NT + nSub × (170N/nSub + 500(N/nSub)^1.5) bytes

## MATLAB Figure Quality

- **Always use LaTeX interpreter** for titles, labels, legends: `'Interpreter', 'latex'`
- Use `$...$` for inline math: `title('$\kappa(M^{-1}A)$ vs $H/\delta$', 'Interpreter', 'latex')`
- Use `\partial`, `\Omega`, `\Gamma`, `\alpha`, `\delta`, `\kappa`, `\rho` etc.
- Set `'Interpreter', 'latex'` on: `title()`, `xlabel()`, `ylabel()`, `legend()`, `text()`
- Use `\setminus` for set difference: `$\partial\Omega_i \setminus \partial\Omega$`

## iFEM Coding Style (Chen Long)

This project follows the sparse matrixlization style from Long Chen's iFEM package:

- **Assemble in one shot** — Build `(ii, jj, ss)` index/value vectors across all elements, then call `sparse(ii, jj, ss, N, N)` once. Never loop over elements assigning into a sparse matrix.
- **Vectorize across elements** — Use element-wise array operations (`.`, `.*`, `./`). Element geometry (area, gradients) is computed once for all elements.
- **Pre-allocate index arrays** — `ii = zeros(nEntries,1); jj = zeros(nEntries,1); ss = zeros(nEntries,1); idx = 0;` — fill in blocks, then truncate.
- **Struct-based API** — Geometry in `mesh` struct (`node`, `elem`, `bdFlag`, `area`, `edge`). PDE data in `pde` struct (`coef`, `source`, `dirichlet`, `neumann`).
- **Function naming** — `assembleXxx` for matrix assembly, `camelCase` for utilities. No underscores in function names (MATLAB convention).
- **One short doc line** — One comment line above the function stating the mathematical formula. No verbose docstrings.
- **Use built-in `sparse` summation** — Duplicate `(i,j)` entries are automatically summed by `sparse`, so you can write the same `(i,j)` pair multiple times and get the accumulated value.

## Reuse-First Principle

Before writing any new function:
1. **Grep** the codebase for existing utilities that already do the job or can be extended.
2. **Check** `verify/` for existing test patterns to reuse.
3. **Use** existing mesh utilities (`edgeMesh2D`, `edgeMesh3D`, `faceMesh3D`, `extendMesh2D`, `extendMesh3D`, `quadtriangle`, `quadtet`).
4. **Use** existing basis evaluators (`lagrange2D`, `lagrange3D`, `nedelec1_2D`, `nedelec1_3D`, `nedelec2_2D`, `nedelec2_3D`).
5. **Prefer extending** an existing function over creating a new one alongside it.
6. **Do not reimplement** quadrature, basis gradients, or mesh topology — they already exist.

## Git Commit Policy

- **Commit when a phase is complete and verified** — After writing a component and its verification passes, commit immediately. Don't batch unrelated changes.
- **Document genuine bugs in commit messages** — When a non-obvious bug was encountered and fixed during development, describe it in the commit body:
  - What the symptom was (wrong output, crash, assertion failure)
  - The root cause (sign error, indexing mistake, missing edge case)
  - How the fix resolves it
  - Format:
    ```
    Fix NE_2: higher-order edge DOF sign parity

    Bug: interior DOF sign was not set to 1, causing sign flips on
    elements with reversed edges via gSign propagation.
    Root cause: gSign(:, 7:8) was left as zeros instead of ones.
    Fix: set gSign(:, 7:8) = 1 after the edge-DOF loop.
    ```
- **Do NOT commit** half-finished work or code that hasn't been verified.

## File Organization

After completing a phase, organize new files into their appropriate folders. Create new folders as needed.

| Folder | Purpose |
|--------|---------|
| `src/Assembly/Lagrange/` | Scalar Lagrange FE assembly routines |
| `src/Assembly/Nedelec/` | Nedelec FE assembly routines |
| `src/FE/Lagrange/` | Lagrange basis and mesh-extension utilities |
| `src/FE/Nedelec/` | Nedelec basis, orientation, and DOF utilities |
| `src/Utils/` | Mesh, quadrature, transfer, and other auxiliary utilities |
| `src/DDM/` | Domain decomposition partitioning and solver routines |
| `src/Preconditioners/` | AS/OAS/ORAS preconditioner builders |
| `verify/` | Numerical verification and test scripts (`verify_*.m`) |
| `debug/` | One-off debugging and investigation scripts (`debug_*.m`) |
| `.Codex/` | Codex configuration (skills, commands, hooks) |

- **Test scripts always go in `verify/`** — e.g., `verify/verify_ned2_2D.m`.
- **Debug/investigation scripts always go in `debug/`** — e.g., `debug/debug_cond.m`.
- **Create a new subfolder under `src/`** when a logical group of library files warrants it.
- **Never leave standalone scripts at root** — they belong in `verify/`, `debug/`, or a topic folder.

## DDM: Overlap Parameter Rule

- **δ must be an integer multiple of h** (δ = k·h, k ∈ ℕ) unless explicitly specified otherwise.
- This ensures subdomain boundaries align with mesh edges, producing straight (non-zig-zag) interfaces on structured meshes.
- Applies to both strip and checkerboard partitions in all dimensions.

## DDM: Mathematical Formulation

### Spaces

| Space | Definition | Boundary Condition |
|-------|-----------|-------------------|
| Fine space V_h | P1 FEM on uniform mesh of size h, dim ≈ h^{-d} | u=0 on ∂Ω |
| Subdomain Ω_i (ASM) | Overlapping: elements with centroid x ∈ [a+(i-1)H-δ, a+iH+δ] | u=0 on ∂Ω_i (Dirichlet inner BC) |
| Subdomain Ω_i^0 (OSM) | Non-overlapping: elements with centroid x ∈ [a+(i-1)H, a+iH) | Robin on Γ_{ij} |
| Interior nodes V_{h,i} | Nodes where ALL incident elements ∈ Ω_i | Free DOFs for subdomain solve |

### ASM (Additive Schwarz) — Overlapping, Dirichlet inner BC

M^{-1} = Σ_i R_i^T A_i^{-1} R_i

- R_i: V_h → V_{h,i} — restricts global free DOFs to interior free DOFs of Ω_i
- A_i = R_i A R_i^T — extracted from global stiffness (no re-assembly needed)
- V_{h,i} = {v ∈ V_h|_{Ω_i} : v = 0 on ∂Ω_i} — Dirichlet on artificial boundaries
- κ(M^{-1}A) ~ O(1 + H/δ), independent of h
- **Overlap is ESSENTIAL**: without overlap (δ=0), Dirichlet inner BC makes subdomains disconnected → κ→∞

### OSM (Optimized Schwarz) — Non-overlapping, Robin transmission

- ∂u_i/∂n_i + α u_i = ∂u_j/∂n_i + α u_j on Γ_{ij}
- Robin term: α·M_Γ added to stiffness, M_Γ·g added to RHS
- Flux from neighbor: g(k) = b_j(k) - (A_j u_j)(k) + α u_j(k)
- Optimal α ≈ 0.5·π/H (empirically)
- ρ independent of h, degrades as H→0: ρ→1 without coarse space

## DDM Verification Commands

### Full parameter study (1D + 2D + 3D)

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_study.m');"
```

Single script covering all dimensions: ASM κ vs h/H/δ, OSM ρ vs α/H, ASM vs OSM comparison.

### Generate partition diagrams

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_partition_viz.m');"
```

Generates in `verify/`:
- `fig_asm_overlap.png` — 2D ASM: elements colored by Ω_i, overlap region, interior★ vs boundary○ nodes
- `fig_osm_nonoverlap.png` — 2D OSM: non-overlapping subdomains, interface edges in red

## Verified Results Summary

### 1D Results (h=1/32..1/256)

| Method | nSub | δ>0? | κ | ρ | Iters |
|--------|------|------|---|---|-------|
| ASM | 2 | yes | 1.0 | ~0 | 3 |
| ASM | 4 | yes | 1.2–1.3 | 0.002–0.004 | 5 |
| ASM | 8 | yes | 2.0–2.5 | 0.030–0.051 | 9 |
| ASM | any | **no** | 80–550 | 1.2–1.5 | **fails** |
| OSM | 2 | — | — | 0.15 | 9 |
| OSM | 4 | — | — | 0.46 | 18 |
| OSM | 8 | — | — | 0.91 | 100 |

### 2D Results (h=1/12..1/24)

| Method | nSub | δ>0? | κ | ρ | Iters |
|--------|------|------|---|---|-------|
| ASM | 2 | yes | 1.4–2.0 | 0.006–0.030 | 6–7 |
| ASM | 3 | yes | 2.6–5.2 | 0.056–0.151 | 9–13 |
| ASM | 4 | yes | 3.8–9.1 | 0.105–0.253 | 11–17 |
| ASM | any | **no** | 700–1600 | 1.1 | **fails** |
| OSM | 2 | — | — | 0.94 | >100 |
| OSM | 3 | — | — | 0.99 | >100 |

### 3D Results (h=1/6..1/8)

| Method | nSub | δ>0? | κ | ρ | Iters |
|--------|------|------|---|---|-------|
| ASM | 2 | yes | 2.0–2.5 | 0.029–0.050 | 7–8 |
| ASM | 3 | yes | 5.3 | 0.155 | 13 |
| ASM | any | **no** | 8000+ | 1.0+ | **fails** |
| OSM | 2 | — | — | 0.96 | >100 |

### Key Findings

1. **ASM-PCG is effective across all dimensions** (3–17 iterations with proper overlap)
2. **Overlap δ > 0 is essential** for ASM with Dirichlet inner BC — without overlap the preconditioned system becomes singular (subdomains disconnected at shared boundaries)
3. **κ is independent of h**, depends on H/δ ratio — confirmed in 1D/2D/3D
4. **OSM converges in 1D** (ρ≈0.15 with optimal α) but **fails in 2D/3D** (ρ→1) without a coarse space
5. **One-level OSM is not viable in 2D/3D** — a two-level method with coarse space is required

### Run DDM extension study (checkerboard, overlapping OSM, coarse space)

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_extension.m');"
```

5-table study on 2D Poisson:
- Table 1: ASM strip vs checkerboard (κ, PCG iters)
- Table 2: OSM strip vs checkerboard (ρ, iters)
- Table 3: Overlapping OSM δ-effect (ρ vs δ, α)
- Table 4: Two-level ASM with coarse space
- Table 5: Two-level OSM with coarse space

**Extension Results Summary:**
- Checkerboard partitioning works with ASM (κ slightly higher than strips for same element count)
- Overlapping OSM: δ=0.16 gives ρ=0.04 (6 iters), vs ρ=0.94 (fails) for non-overlapping. Overlap dramatically improves OSM.
- Two-level OSM: P1 coarse H=1/6 gives ρ=0.17 (9 iters), vs ρ=0.996 (fails) for one-level. Coarse space transforms OSM from useless to practical in 2D.
- Two-level ASM: additive coarse correction increases κ when fine preconditioner is already strong. Multiplicative/hybrid needed for benefit.
- P2 DDM: requires special DOF classification (vertex vs edge interior/boundary) — open problem.

## ORAS for Helmholtz (Gong-Gander-Graham-Spence 2021/2022)

### Mathematical Formulation

Helmholtz with impedance BC: `-(Δ + k²)u = f`, `∂u/∂n - iku = g` on ∂Ω.

Sesquilinear form: `a(u,v) = ∫_Ω(∇u·∇v̄ - k²uv̄)dx - ik∫_{∂Ω}uv̄ ds`

ORAS preconditioner: `B_h^{-1} = Σ_j R̃_{h,j} A_{h,j}^{-1} R_{h,j}`

- `R_{h,j}`: restriction (by duality)
- `R̃_{h,j}`: weighted prolongation = `R_{h,j}(I_h(χ_j · v))` — multiply local solution by partition of unity χ_j before extending
- `A_{h,j}`: local Helmholtz with impedance BC on ALL of ∂Ω_j
- χ_j(node) = 1/(#subdomains containing node)

GMRES is the preferred global solver (non-Hermitian A). Richardson needs strong damping.

### Usage: ORAS with adjustable k, h, delta, and solver mode

```matlab
% ---- Parameters -----------------------------------------------------------
k     = 20;                   % wavenumber
h     = 2/(10*k);             % mesh size (h*k = 0.2 for P1 resolution)
delta = 1/4;                  % overlap extension per side (physical distance)
nSub  = 8;                    % number of strip subdomains
% For checkerboard: nSub = [nx, ny], e.g. [5, 5]

% ---- Mesh -----------------------------------------------------------------
[node, elem, bd] = squaremesh([0, 16/3, 0, 1], h);  % strip domain
% [node, elem, bd] = squaremesh([0, 1, 0, 1], h);    % checkerboard domain

% ---- Manufactured solution: u = sin(pi x) sin(pi y) -----------------------
u_ex = @(x,y) sin(pi*x) .* sin(pi*y);
f    = @(x,y) (2*pi^2 - k^2) * u_ex(x,y);
g    = @(x,y) 0;              % impedance BC: du/dn - iku = g

% ---- Global Helmholtz assembly -------------------------------------------
[A, b] = assembleHelmholtz2D(node, elem, bd, k, f, g);

% ---- Partition ------------------------------------------------------------
parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);

% ---- Build ORAS preconditioner -------------------------------------------
% solverMode: 'lu' (default, fast iterations, more memory — may exceed RAM)
%             'direct' — A\b each iteration via UMFPACK, less memory
applyPrecon = orasHelmholtz(node, elem, bd, k, parts, 1, 'direct');

% ---- Solve with Richardson iteration -------------------------------------
u = zeros(size(A,1), 1);
for it = 1:200
    r = b - A * u;
    u = u + applyPrecon(r);
    if norm(r)/norm(b) < 1e-6, break; end
end
% ---- Or solve with GMRES -------------------------------------------------
% [u, flag, relres, iter] = gmres(A, b, [], 1e-6, 200, applyPrecon);
```

**Parameter guide:**
| Param | Description | Typical value |
|-------|-------------|---------------|
| `k` | Wavenumber | 20 (h=0.01), 40 (h=0.005), etc. |
| `h` | Mesh size | `2/(10*k)` for h·k=0.2 resolution |
| `delta` | Overlap per side | `1/4` for strips, `H/4` for checkerboard |
| `nSub` | Subdomain count | Scalar→strips, `[nx,ny]`→checkerboard |
| `degree` | FE order | 1 (P1), 2 (P2 experimental) |
| `solverMode` | Subdomain solver | `'lu'` (fast, high memory), `'direct'` (backslash, low memory) |
| `useParfor` | Parallel subdomain setup | `false` (default), `true` (needs active `parpool`) |

**Solver mode comparison:**
| Mode | Peak memory | Per-iteration speed | When to use |
|------|------------|---------------------|-------------|
| `'lu'` | High (LU fill-in can be 10-100× matrix size) | Fast (forward/back substitution) | Small/medium problems, enough RAM |
| `'direct'` | Low (stores only A_j, no factors) | Slower (UMFPACK each iteration) | Large 2D/3D problems, HPC |

### Run overlapping OSM comprehensive study

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ddm_osm_overlap.m');"
```

Shows partition diagrams first, then 4 tables:
- Table 1: Overlapping OSM — δ effect for strip partitions (nSub=2,3,4)
- Table 2: Overlapping OSM — δ effect for checkerboard (2×2, 3×3)
- Table 3: Overlapping OSM — α sensitivity scan (best δ=H/2)
- Table 4: Two-level overlapping OSM — coarse + overlap combined (strip + checkerboard)

Generates diagrams:
- `fig_osm_strip_overlap.png` — strip non-overlap vs overlap comparison
- `fig_osm_checkerboard_overlap.png` — checkerboard 2×2 and 3×3 with overlap

### Show subdomain boundaries (zig-zag visualization)

```bash
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/plot_subdomain_boundaries.m');"
```

Generates `verify/fig_smooth_boundaries.png` — 3 panels contrasting zig-zag vs straight boundaries:
- δ not aligned with mesh → zig-zag
- δ = h → straight vertical boundaries
- δ = 2h → straight, wider overlap
- **Rule:** choose δ = k·h (integer multiple of mesh size) for straight subdomain boundaries on structured meshes.
