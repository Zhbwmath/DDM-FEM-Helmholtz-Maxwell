# DDM-FEM-Helmholtz-Maxwell Project

Created: 2026-05-21
Updated: 2026-05-25

## MATLAB Execution

- **Always run MATLAB silently**: use `-nosplash -nodesktop -batch` — no console windows.
- **Never use `-noFigureWindows`** unless the user explicitly asks to suppress all graphics. Experiment figures must display.
- For running a script: `matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('script.m');"`
- For running inline code: `matlab -nosplash -nodesktop -batch "addpath(genpath('.')); <code>;"`
- Check for MATLAB at `C:\Program Files\MATLAB\R2023a\bin\matlab.exe` first, then fall back to `matlab` on PATH.

## HPC Rules (Workstation: 48 cores, 549 GB RAM)

- **Ask permission before any test estimated to use >200 GB memory.** Estimate memory usage and present it before running.
- **Use `parfor` for subdomain setup.** Start parpool before large runs: `parpool('local', feature('numcores'))`.
- **Vectorize all assembly.** Use iFEM-style one-shot `sparse(ii, jj, ss, N, N)`.
- **Avoid `containers.Map` with string keys for large meshes.** Use `ismember` on sorted edge lists or sparse hashing.
- **Memory estimation rule of thumb** for 2D Helmholtz with N nodes, NT elements, nSub subdomains:
  - Global sparse matrix A: ~7N nonzeros × 16 bytes (complex) ≈ 112N bytes
  - Edge list: ~6NT entries × 8 bytes ≈ 48NT bytes
  - Per subdomain (N/nSub nodes): ~3 × 7×N/nSub nonzeros × 8 bytes (real) + LU ≈ 500×(N/nSub)^1.5 bytes
  - Total ≈ 112N + 48NT + nSub × (170N/nSub + 500(N/nSub)^1.5) bytes

## Wave-Problem Resolution Rules

- **Fine mesh for Helmholtz/Maxwell verification must resolve pre-asymptotic accuracy.** In addition to any coarse-mesh condition such as `H*k = O(1)`, the fine/reference Lagrange mesh for wave problems must satisfy at least `h = O(k^(-(2*p+1)/(2*p)))`, where `p` is the polynomial degree. For P1 this is `h = O(k^(-3/2))`. LOD experiments must choose the reference/fine mesh from this rule before judging pollution or localization behavior.
- **Record the fine-scale rule in reproduction notes.** Any Helmholtz/Maxwell LOD or DDM reproduction report must state the selected `p`, `k_max`, resulting `h`, and whether the fine mesh is nested with every coarse mesh in the sweep.
- **Modify experiment records in place.** When follow-up requirements change an existing experiment, update the previous reproduction `.md` file directly instead of creating a new report. If a table already records the affected data, modify the noncompliant columns or rows in that table rather than adding a replacement table. This keeps experiment documentation minimal and preserves continuity.

## MATLAB Figure Quality

- **Always use LaTeX interpreter** for titles, labels, legends: `'Interpreter', 'latex'`
- Use `$...$` for inline math: `title('$\kappa(M^{-1}A)$ vs $H/\delta$', 'Interpreter', 'latex')`
- Use `\partial`, `\Omega`, `\Gamma`, `\alpha`, `\delta`, `\kappa`, `\rho` etc.
- Set `'Interpreter', 'latex'` on: `title()`, `xlabel()`, `ylabel()`, `legend()`, `text()`
- Use `\setminus` for set difference: `$\partial\Omega_i \setminus \partial\Omega$`

## Documentation Metadata

- Every new or substantially updated Markdown document must include `Created: YYYY-MM-DD` and `Updated: YYYY-MM-DD` near the top. Keep `Created` fixed and refresh `Updated` when the document is changed.
- Reproduction documents must also include `Verification entry point:` with the rerunnable script/function/command, and `Main utilities:` listing the principal assembly, solver, preconditioner, mesh, or verification functions used.
- Active research notes under `tasks/<topic>/` follow the same metadata rule. When a task note is promoted to `docs/`, preserve the original creation date, refresh the update date, and keep the verification entry point current.

## Paper Resources

- Store user-provided papers for local reading under `resources/`.
- Do not commit paper PDFs or other copyrighted source documents. Commit only lightweight metadata, notes, or directory marker files that explain what should live there.

## iFEM Coding Style (Chen Long)

This project follows the sparse matrixlization style from Long Chen's iFEM package:

- **Assemble in one shot** — Build `(ii, jj, ss)` index/value vectors across all elements, then call `sparse(ii, jj, ss, N, N)` once. Never loop over elements assigning into a sparse matrix.
- **Vectorize across elements** — Use element-wise array operations (`.`, `.*`, `./`). Element geometry (area, gradients) is computed once for all elements.
- **Vectorize edge/face terms too** — For jump, trace, boundary, and interface integrals, collect all relevant edges/faces first and evaluate geometry, orientations, quadrature traces, and jump values as arrays. Avoid loops over edges/faces in production assemblers; small loops over fixed quadrature points, polynomial degree, or derivative order are acceptable when the edge/element dimension is vectorized.
- **Keep loop references for refactors** — When replacing a clear edge/element loop prototype with a vectorized assembler, add a focused verification that compares the new sparse matrix against the loop reference before relying on stronger mathematical tests.
- **Pre-allocate index arrays** — `ii = zeros(nEntries,1); jj = zeros(nEntries,1); ss = zeros(nEntries,1); idx = 0;` — fill in blocks, then truncate.
- **Struct-based API** — Geometry in `mesh` struct (`node`, `elem`, `bdFlag`, `area`, `edge`). PDE data in `pde` struct (`coef`, `source`, `dirichlet`, `neumann`).
- **Function naming** — `assembleXxx` for matrix assembly, `camelCase` for utilities. No underscores in function names (MATLAB convention).
- **One short doc line** — One comment line above the function stating the mathematical formula. No verbose docstrings.
- **Use built-in `sparse` summation** — Duplicate `(i,j)` entries are automatically summed by `sparse`, so you can write the same `(i,j)` pair multiple times and get the accumulated value.

## Numerical Coding Discipline

These rules adapt the coding principles in `Karpathy's CLAUDE.md` to this FEM/DDM MATLAB project.

### Think Before Coding

- State the mathematical interpretation before implementation: PDE, weak form, finite element space, boundary conditions, and matrix/operator form.
- If a paper, existing code path, or user instruction is ambiguous, ask before coding. Do not fill gaps with plausible numerical-analysis folklore.
- Surface important tradeoffs early: paper fidelity vs. memory, direct LU vs. iterative solve, exact table parameters vs. scaled verification.
- For DDM work, clarify the geometry and interface conditions before writing partition or solver code.

### Simplicity First

- Implement the smallest paper-faithful component that can be verified.
- Do not add general frameworks, unused options, or speculative solver modes unless required by the paper or requested by the user.
- Prefer extending existing assembly, mesh, quadrature, partition, and preconditioner utilities over adding parallel versions.
- If a verification script becomes broad or slow, split focused checks from paper-scale reproduction instead of mixing them.

### Surgical Changes

- Touch only the files needed for the requested method, reproduction, or verification.
- Preserve existing MATLAB style: vectorized iFEM-style assembly, one-shot `sparse(ii,jj,ss,...)`, camelCase utilities, short formula comments.
- Do not refactor unrelated code while implementing a paper method. Mention unrelated issues separately.
- Remove only unused code introduced by the current change; do not clean up pre-existing dead code unless asked.

### Goal-Driven Verification

For each implementation/reproduction phase, define success criteria before running:

1. Formulation extraction -> verify: equations and algorithm steps cite the paper section/equation/table.
2. Matrix translation -> verify: dimensions, DOF sets, restrictions, local operators, and boundary terms are explicit.
3. Implementation -> verify: focused numerical checks pass on small meshes.
4. Paper reproduction -> verify: generated table is compared directly against the paper table.
5. Closeout -> verify: scripts live in `verify/`, debug helpers in `debug/`, and deviations from the paper are documented.

A task is complete only when the result table answers the reproduction question: consistent with the paper, inconsistent, or blocked by a clearly stated limitation.

## Reuse-First Principle

Before writing any new function:
1. **Grep** the codebase for existing utilities that already do the job or can be extended.
2. **Check** `verify/` for existing test patterns to reuse.
3. **Use** existing mesh utilities (`edgeMesh2D`, `edgeMesh3D`, `faceMesh3D`, `extendMesh2D`, `extendMesh3D`, `quadtriangle`, `quadtet`).
4. **Use** existing basis evaluators (`lagrange2D`, `lagrange3D`, `nedelec1_2D`, `nedelec1_3D`, `nedelec2_2D`, `nedelec2_3D`).
5. **Prefer extending** an existing function over creating a new one alongside it.
6. **If no utility exists, create one for reuse** — Do not bury reusable numerical pieces such as normal-derivative jumps, trace matrices, restrictions, or edge/face geometry inside one paper-specific assembler. Extract a small reusable subroutine and have the specific assembler call it.
7. **Do not reimplement** quadrature, basis gradients, or mesh topology — they already exist.

## Research Subagents

- **Use `math-searcher`** (`.claude/agents/math-searcher.md`) when a request needs internet literature search, article extraction, or implementation search for DDM, FEM, Helmholtz, or Maxwell topics.
- Give `math-searcher` a bounded target: method names, equations/sections to extract, desired source type (paper, arXiv, code, documentation), and implementation language if relevant.
- `math-searcher` should prioritize primary sources, return URLs/DOIs/arXiv IDs, extract only the requested formulas or algorithm details, and state how each result maps to this MATLAB codebase.
- Do not treat internet summaries as implementation authority. Convert any extracted formulation into this project's notation and verify locally before coding.
- **Use `math-translator`** for paper reproduction and active research tasks after the source formulation is identified. It should write the PDE, boundary/interface conditions, variational form, integration-by-parts steps when they matter, and matrix/operator representation into the task's Markdown file under `tasks/<topic>/` or the relevant `docs/<article-or-topic>/` folder.
- When several mathematically equivalent discretizations exist, `math-translator` must name the alternatives and state which one this repo uses for the task. For example, PML notes should explicitly say whether the implementation uses a divergence-form stretched-coordinate bilinear form or an expanded non-divergence form, then give the corresponding matrix formula.
- The task Markdown should be the durable formulation record for ongoing research work: update it as implementation choices change instead of leaving the formulas only in chat or temporary scratch notes.

## Paper Reproduction Workflow

Use this workflow whenever the user asks to **reproduce**, **replicate**, or **match** experiments from a paper.

- **Document metadata required:** the first block of each reproduction Markdown file must state the reproduction target, `Created`, `Updated`, `Verification entry point`, and `Main utilities`.
- **Goal first:** the goal is not to improve the method, tune aggressively, or make a new benchmark. The goal is to determine whether this repo can produce tables/figures consistent with the paper.
- **Extract before coding:** use `math-searcher` when needed to find the paper, preprint, author implementation, supplementary material, or related code. Extract the exact algorithm, PDE, boundary conditions, discretization, stopping rules, reported metrics, and table/figure parameters.
- **Translate to matrices:** use `math-translator` when needed to convert the paper formulation into this repo's matrix notation: global operator, local subdomain operators, restriction/prolongation, partition of unity, coarse space, transmission terms, and solver iteration.
- **Parameter sheet required:** write the paper parameters into a concrete experiment form before running:
  - domain and boundary conditions
  - PDE coefficients, wavenumber/frequency, material parameters
  - mesh size `h`, subdomain size `H`, overlap `delta`, polynomial degree, quadrature
  - number/type of subdomains, partition geometry, coarse space
  - solver, preconditioner, tolerance, max iterations, restart/damping
  - reported paper table/figure target values
- **Strict alignment:** match the paper unless a project rule or HPC rule prevents it. Do not silently replace algorithms, boundary conditions, solvers, partitioning, or parameters with convenient alternatives.
- **HPC exception:** if exact paper parameters are estimated to exceed the active HPC permission threshold or are otherwise unsafe, stop and report the memory estimate. Propose the closest scaled experiment separately and label it as scaled, not reproduced.
- **Comparison report:** every reproduction run should end with a compact table comparing `paper value`, `repo value`, `relative/absolute difference`, and `notes`. State whether the result is consistent, partially consistent, or inconsistent.
- **No hidden tuning:** if extra tuning is needed to match the paper, document it as a deviation. Keep the paper-faithful run as the baseline.
- **Temporary reproduction notes:** during an active literature-reproduction run, temporary Markdown notes and generated figures may live under a dedicated `verify/<paper-or-method>/` folder.
- **Article folder naming:** completed article-reproduction folders under `docs/` must use the pattern `<AMS-style citation abbreviation>_<brief method>`, e.g. `GGGLS24_pml`, `TW05_asm_poisson`, or `Gan06_osm_poisson`.
- **One article per reproduction doc:** do not mix reproduction experiments from different articles in one report. Split mixed reports by paper/book/article target. The first line of each reproduction Markdown file must state `Reproduction target: ...`.
- **Interest-driven exception:** exploratory experiments driven by project interests rather than a specific paper do not need the article-abbreviation naming rule, but should still use clear folder names.
- **Closeout cleanup:** before committing a completed reproduction, remove half-finished scratch files from `verify/` and move the finished article-level report folder to `docs/<AMS-style abbreviation>_<brief method>/`. Keep `verify/` for executable checks, temporary run artifacts, and scripts that can be rerun.

## Git Commit Policy

- **Include user document edits when committing** — If the user has modified Markdown, task notes, reports, or project-rule documents related to the current phase, include those document changes in the commit and follow the updated rules without asking again.
- **Double-check user code edits before committing** — If the user has modified source or verification code during the phase, inspect the diff and ask only when the intent is unclear, the code conflicts with the current implementation, or committing it would mix unrelated work.

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
| `tasks/` | Active research-task folders based on the current repo; keep task-local formulation notes, matrix translations, implementation plans, open questions, and intermediate Markdown records here until they become stable documentation |
| `docs/` | Stable documentation, article reproduction reports, result summaries, and generated figures referenced by reports |
| `verify/` | Numerical verification and test scripts (`verify_*.m`) |
| `debug/` | One-off debugging and investigation scripts (`debug_*.m`) |
| `.claude/agents/` | Project sub-agent definitions such as `math-searcher` |
| `.claude/skills/` | Project Claude skills and migrated command helpers |
| `.agents/skills/` | Project Codex skills and source-command wrappers |
| `.Codex/` | Codex configuration (skills, commands, hooks) |

- **Test scripts always go in `verify/`** — e.g., `verify/verify_ned2_2D.m`.
- **Debug/investigation scripts always go in `debug/`** — e.g., `debug/debug_cond.m`.
- **Active research notes go in `tasks/<topic>/`** until they are finished enough to move into `docs/`; keep task-local formulas, matrix translations, parameter sheets, and unresolved implementation choices there.
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

Stable DDM result summaries and ORAS/Helmholtz reproduction notes live in:
- `docs/TW05_asm_poisson/`
- `docs/Gan06_osm_poisson/`
- `docs/GGGS21_oras_helmholtz/`
- `docs/GGGLS24_pml/`
- `docs/GGGLS24_ras_pml/`
- `docs/interest_theory_checks/`
