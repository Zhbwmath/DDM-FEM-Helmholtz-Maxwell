Created: 2026-05-26
Updated: 2026-05-27
Verification entry point: `verify/verify_lod_patch_struct.m`
Main utilities: `buildLOD`, `lodBuildPatches`, `helmholtzLODProblem2D`, `helmholtzLODProblem3D`, `weightedClementP1`

# LOD API Plan

The generic entry point is
```matlab
lod = buildLOD(nodeH, elemH, nodeh, elemh, problem, opts);
```

The `problem` struct supplies matrix callbacks:
```matlab
problem.form.global()                    % returns [A,b]
problem.form.patch(patch,T)              % returns A_P
problem.form.elementRhs(T,dofs,patch,T,P)
problem.form.elementRhsAdjoint(T,dofs,patch,T,P)
problem.constraints.patch(Q,patch,T,opts)
problem.transfer()                       % optional, default prolongateNestedP1
problem.interpolation()                  % optional, default weightedClementP1
```

The matrix convention is always `v' * A * u = a(u,v)`. The constraint callback returns `C_P`, and the production corrector solve enforces `C_P' * q = 0` through the saddle system.

The return value groups data as:
```matlab
lod.basis.coarse
lod.basis.trial
lod.basis.test
lod.system.A
lod.system.b
lod.system.AH
lod.system.bH
lod.patch
lod.solution
```

Patch metadata is aggregated under `lod.patch` with cell fields for `coarseElemIds`, `fineElemIds`, `local2global`, `freeLocalDof`, `boundaryLocalDof`, `artificialBoundaryLocalDof`, `physicalBoundaryLocalDof`, and a `stats` struct array after corrector construction.

Patch construction dispatches by mesh dimension. Triangular P1 meshes use the original 2D path, and tetrahedral P1 meshes use the 3D path with physical boundary DOFs taken from tetrahedral face flags and artificial boundary DOFs fixed out of the local saddle solve.

`opts.useParfor=true` enables basis-wise patch corrector parallelism. The numerical path is unchanged: each worker solves the same local constrained saddle systems.
