% VERIFY_LOD_HELMHOLTZ3D_SMOKE  Small Petrov-Galerkin LOD Helmholtz solve in 3D.

fprintf('========== LOD Helmholtz 3D Smoke Verification ==========\n\n');

[nodeH, elemH, bdH] = cubemesh([0, 1, 0, 1, 0, 1], 1);
[nodeh, elemh, bdh] = cubemesh([0, 1, 0, 1, 0, 1], 0.5);
opts = struct('oversampling', 1, 'solveCoarse', true);
lod = buildLODHelmholtz3D(nodeH, elemH, bdH, nodeh, elemh, bdh, 1, 1, 0, opts);

fprintf('Test 1: 3D coarse Petrov-Galerkin system is finite ... ');
assert(all(isfinite(nonzeros(lod.system.AH))), '3D LOD coarse matrix has non-finite entries.');
assert(all(isfinite(lod.system.bH)), '3D LOD coarse RHS has non-finite entries.');
assert(rcond(full(lod.system.AH)) > 1e-14, '3D LOD coarse matrix is numerically singular.');
fprintf('passed\n');

fprintf('Test 2: 3D LOD solution and fine comparison are finite ... ');
assert(all(isfinite(lod.solution.uh)), '3D LOD solution has non-finite entries.');
assert(isfield(lod.solution, 'relativeEnergyError'), 'Missing 3D relative energy error.');
assert(isfinite(lod.solution.relativeEnergyError), '3D relative energy error is non-finite.');
assert(lod.solution.relativeEnergyError < 2, '3D smoke-scale LOD error is unexpectedly large.');
fprintf('passed (relative energy error %.3e)\n', lod.solution.relativeEnergyError);

fprintf('Test 3: 3D patch stats are populated ... ');
assert(numel(lod.patch.stats) == size(elemH, 1), 'Missing 3D patch stats.');
assert(all([lod.patch.stats.freeDof] > 0), 'Each 3D patch must have free DOFs.');
assert(all(isfinite([lod.patch.stats.primalResidual])), '3D primal residuals are non-finite.');
assert(all(isfinite([lod.patch.stats.adjointResidual])), '3D adjoint residuals are non-finite.');
fprintf('passed\n');

fprintf('\n========== LOD Helmholtz 3D smoke tests PASSED ==========\n');
