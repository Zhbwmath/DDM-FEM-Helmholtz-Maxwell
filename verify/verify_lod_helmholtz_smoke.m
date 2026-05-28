% VERIFY_LOD_HELMHOLTZ_SMOKE  Small Petrov-Galerkin LOD Helmholtz solve.

fprintf('========== LOD Helmholtz Smoke Verification ==========\n\n');

[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);
opts = struct('oversampling', 1, 'solveCoarse', true);
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, nodeh, elemh, bdh, 3, 1, 0, opts);

fprintf('Test 1: coarse Petrov-Galerkin system is finite ... ');
assert(all(isfinite(nonzeros(lod.system.AH))), 'LOD coarse matrix has non-finite entries.');
assert(all(isfinite(lod.system.bH)), 'LOD coarse RHS has non-finite entries.');
assert(rcond(full(lod.system.AH)) > 1e-14, 'LOD coarse matrix is numerically singular.');
fprintf('passed\n');

fprintf('Test 2: LOD solution and fine comparison are finite ... ');
assert(all(isfinite(lod.solution.uh)), 'LOD solution has non-finite entries.');
assert(isfield(lod.solution, 'relativeEnergyError'), 'Missing relative energy error.');
assert(isfinite(lod.solution.relativeEnergyError), 'Relative energy error is non-finite.');
assert(lod.solution.relativeEnergyError < 2, 'Smoke-scale LOD error is unexpectedly large.');
fprintf('passed (relative energy error %.3e)\n', lod.solution.relativeEnergyError);

fprintf('\n========== LOD Helmholtz smoke tests PASSED ==========\n');
