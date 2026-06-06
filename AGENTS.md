# DDM-FEM-Helmholtz-Maxwell Project

Created: 2026-05-21
Updated: 2026-06-07

## MATLAB Execution

- **Prefer the MATLAB MCP server when it is installed and healthy.** For MATLAB inspection, small code evaluation, scripts, tests, and Code Analyzer checks, use the available MCP tools first (`detect_matlab_toolboxes`, `evaluate_matlab_code`, `run_matlab_file`, `run_matlab_test_file`, `check_matlab_code`). Run from the repo root or pass the project path so the session executes with `addpath(genpath('.'))` when repo utilities are needed.
- **Smoke-test MCP before relying on it.** A successful `detect_matlab_toolboxes` or a small `evaluate_matlab_code` call is enough to treat MCP as working for the current session. If the MCP server is missing, unresponsive, using the wrong MATLAB installation, or fails before MATLAB code starts, record that failure and use the batch fallback below.
- **Batch fallback:** run MATLAB silently with `-nosplash -nodesktop -batch` and no console windows. For a script, use `matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('script.m');"`. For inline code, use `matlab -nosplash -nodesktop -batch "addpath(genpath('.')); <code>;"`.
- **Never use `-noFigureWindows`** unless the user explicitly asks to suppress all graphics. Experiment figures must display.
- For batch fallback, check for MATLAB at `C:\Program Files\MATLAB\R2023a\bin\matlab.exe` first, then fall back to `matlab` on PATH.
- **License errors in sandbox:** if MCP or batch MATLAB fails with a license checkout/license manager error from the Codex sandbox, treat it as an environment/licensing problem, not a numerical-code failure. Report the exact license error text or license manager code, then retry the same MATLAB action outside the sandbox only with the required user approval. Do not edit license files, hard-code license servers, expose license numbers in committed files, or change repo code to work around licensing. If approval is denied or the license still fails outside the sandbox, mark MATLAB verification as blocked by licensing and continue only with non-executing checks.

## HPC Rules (Workstation: 48 cores, 549 GB RAM)

- **Ask permission before any test estimated to use >200 GB memory.** Estimate memory usage and present it before running.
- **Use `parfor` for subdomain setup.** Start parpool before large runs: `parpool('local', feature('numcores'))`.
- **Vectorize all assembly.** Use iFEM-style one-shot `sparse(ii, jj, ss, N, N)`.
- **Memory estimation rule of thumb** for 2D Helmholtz with `N` nodes, `NT` elements, `nSub` subdomains, and GMRES restart/basis length `m`:
  - Global sparse matrix `A`: about `7N` nonzeros x 16 bytes (complex) = `112N` bytes.
  - Edge/topology work arrays: about `6NT` integer entries x 8 bytes = `48NT` bytes.
  - Local sparse matrices: for local size `n_l`, store stiffness/mass/boundary matrices and assembly work as about `3 x 7n_l` complex nonzeros = `336n_l` bytes.
  - Sparse LU factor storage: do **not** estimate as dense. For 2D sparse direct solves, use `16*c_lu*n_l*log2(max(n_l,2))` bytes with default fill constant `c_lu=20`; use `c_lu=30--40` for a conservative indefinite Helmholtz estimate or when ordering/fill behavior is unknown.
  - GMRES Arnoldi basis: about `16N(m+3)` bytes for complex vectors, plus small Hessenberg/orthogonalization work arrays.
  - One-level DDM sparse-LU total:
    `112N + 48NT + nSub*(336n_l + 16*c_lu*n_l*log2(max(n_l,2))) + 16N(m+3)` bytes, where `n_l` should be the actual patch node count if known; otherwise use `N/nSub` only as a coarse first estimate.
  - MATLAB direct/backslash fallback: when using direct solves without storing reusable local LU factors, use the sparse-LU estimate as a peak-memory guide but allow a looser safety margin because MATLAB may release temporaries between solves. Log observed peak memory when possible.
  - Dense fallback estimate, only for a worst-case warning, is `16*n_l^2` bytes per local dense factor. Do not use the old `500*n_l^1.5` term as a memory estimate; it is closer to sparse-direct work scaling than storage and can greatly overestimate memory.

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
- **Markdown math rendering:** use `$...$` for inline formulas and `$$...$$` for display formulas in Markdown files and chat derivations. Avoid `\(...\)` and `\[...\]` because they may not render consistently in this workspace.


## iFEM Coding Style 

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
- **Do not commit**  paper PDFs or other copyrighted source documents. Commit only lightweight metadata, notes, or directory marker files that explain what should live there.

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
| `resources` | Store user-provided papers for local reading |
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

- **DDM mesh-alignment principle:** `delta`, each overlapping subdomain size `H_l` for `Omega_l`, and any two-level coarse mesh size `H` must be integer multiples of the fine mesh size `h` unless the user explicitly specifies an exception.
- In particular, use `delta = k_delta h`, `H_l = k_l h`, and, when a coarse space is used, `H = k_H h` with integer `k_delta`, `k_l`, and `k_H`.
- This ensures overlap boundaries, subdomain boundaries, and coarse-space mesh boundaries align with the fine mesh, producing straight (non-zig-zag) interfaces on structured meshes.
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
- A_i = R_i A R_i^T — extracted from global stiffness
- V_{h,i} = {v ∈ V_h|_{Ω_i} : v = 0 on ∂Ω_i} — Dirichlet on artificial boundaries
- **Overlap is ESSENTIAL**: without overlap (δ=0), Dirichlet inner BC makes subdomains disconnected → κ→∞

### OSM (Optimized Schwarz) — Non-overlapping, Robin transmission

- ∂u_i/∂n_i + α u_i = ∂u_j/∂n_i + α u_j on Γ_{ij}
- Flux from neighbor
- ρ independent of h, degrades as H→0: ρ→1 without coarse space

## DDM Verification Commands

Stable DDM result summaries and ORAS/Helmholtz reproduction notes live in:
- `docs/TW05_asm_poisson/`
- `docs/Gan06_osm_poisson/`
- `docs/GGGS21_oras_helmholtz/`
- `docs/GGGLS24_ras_rms_pml/`

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **DDM-FEM-Helmholtz-Maxwell** (522 symbols, 513 relationships, 0 execution flows as of 2026-06-05). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

> If any GitNexus tool warns the index is stale, refresh the snapshot before relying on graph results. On this Windows/Codex workstation, do **not** assume `npx`, `node` on `PATH`, or `gitnexus.cmd` works: `npx` may be absent, and `node` on `PATH` can resolve to the Codex WindowsApps executable and fail with `Access is denied`.
>
> Use the verified direct invocation from the repo root:
> `C:\Users\Administrator\Documents\Codex\tools\node-v24.14.0-win-x64\node.exe C:\Users\Administrator\Documents\Codex\tools\gitnexus-cli-npm\node_modules\gitnexus\dist\cli\index.js analyze --index-only .`
>
> If the refresh fails with `EPERM` on `C:\Users\Administrator\.gitnexus\registry.json`, rerun the same direct invocation with the required filesystem approval because GitNexus must update its global registry outside the workspace. The local Codex hook path config should live at `C:\Users\Administrator\.codex\hooks\.local\gitnexus.paths.ps1` and point to the same Node and GitNexus CLI files.

## Always Do

- **MUST run impact analysis before editing any symbol.** Before modifying a function, class, or method, run `gitnexus_impact({target: "symbolName", direction: "upstream"})` and report the blast radius (direct callers, affected processes, risk level) to the user.
- **MUST run `gitnexus_detect_changes()` before committing** to verify your changes only affect expected symbols and execution flows.
- **MUST warn the user** if impact analysis returns HIGH or CRITICAL risk before proceeding with edits.
- When exploring unfamiliar code, use `gitnexus_query({query: "concept"})` to find execution flows instead of grepping. It returns process-grouped results ranked by relevance.
- When you need full context on a specific symbol — callers, callees, which execution flows it participates in — use `gitnexus_context({name: "symbolName"})`.

## Never Do

- NEVER edit a function, class, or method without first running `gitnexus_impact` on it.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis.
- NEVER rename symbols with find-and-replace — use `gitnexus_rename` which understands the call graph.
- NEVER commit changes without running `gitnexus_detect_changes()` to check affected scope.

## Resources

| Resource | Use for |
|----------|---------|
| `gitnexus://repo/DDM-FEM-Helmholtz-Maxwell/context` | Codebase overview, check index freshness |
| `gitnexus://repo/DDM-FEM-Helmholtz-Maxwell/clusters` | All functional areas |
| `gitnexus://repo/DDM-FEM-Helmholtz-Maxwell/processes` | All execution flows |
| `gitnexus://repo/DDM-FEM-Helmholtz-Maxwell/process/{name}` | Step-by-step execution trace |

## CLI

| Task | Read this skill file |
|------|---------------------|
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus/gitnexus-cli/SKILL.md` |

<!-- gitnexus:end -->
