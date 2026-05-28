% VERIFY_LOD_CORRECTOR_ORTHOGONALITY  Check local LOD saddle residuals.

fprintf('========== LOD Corrector Orthogonality Verification ==========\n\n');

[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);
opts = struct('oversampling', 1, 'solveCoarse', true);
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, nodeh, elemh, bdh, 2, 1, 0, opts);

fprintf('Test 1: primal and adjoint saddle residuals ... ');
pr = [lod.patch.stats.primalResidual];
ar = [lod.patch.stats.adjointResidual];
assert(all(pr < 1e-9), 'Primal saddle residual too large.');
assert(all(ar < 1e-9), 'Adjoint saddle residual too large.');
fprintf('passed\n');

fprintf('Test 2: primal and adjoint constraints ... ');
cr = [lod.patch.stats.constraintResidual];
acr = [lod.patch.stats.adjointConstraintResidual];
assert(all(cr < 1e-9), 'Primal constraint residual too large.');
assert(all(acr < 1e-9), 'Adjoint constraint residual too large.');
fprintf('passed\n');

fprintf('\n========== LOD corrector orthogonality tests PASSED ==========\n');
