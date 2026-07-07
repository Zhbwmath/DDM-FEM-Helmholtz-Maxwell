% VERIFY_PML_LXZZ_HYBRID_INSTANCE  Complete PML-LOD/LXZZ local-PML smoke.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

fprintf('========== PML LOD LXZZ hybrid instance check ==========\n');

k = 3;
pmlBox = [0, 1, 0, 1];
physicalBox = [0.25, 0.75, 0.25, 0.75];
[node, elem, bdFlag] = squaremesh(pmlBox, 0.25);
[nodeH, elemH, bdH] = squaremesh(pmlBox, 0.5);
pml = struct('physicalBox', physicalBox, 'pmlBox', pmlBox, ...
    'sigmaMax', k, 'sigmaOrder', 2, 'quadOrder', 4);

parts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', 0.25);
parts = smoothPartitionOfUnity2D(parts, pmlBox, [2, 2], 0.25);

opts = struct();
opts.lodOptions = struct('oversampling', 1, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', false, ...
    'constraintTolerance', 1e-12);
opts.localOptions = struct('solverMode', 'direct', ...
    'localStorage', 'matrix', 'localPMLMode', 'subdomain', ...
    'applyMode', 'auto');
opts.adjointType = 'reference';

pre = buildPMLLxzzHybridHelmholtz2D(node, elem, bdFlag, k, pml, ...
    parts, nodeH, elemH, bdH, opts);

assert(strcmp(pre.variant, 'pml'), 'Hybrid wrapper did not report PML variant.');
assert(strcmp(pre.localSolver.info.variant, 'pml'), ...
    'Local solver did not report PML variant.');
assert(any(strcmp(pre.localSolver.info.localPMLModes, 'subdomain')), ...
    'Local PML solver did not use subdomain PML boxes.');
assert(size(pre.A, 1) == numel(pre.fineSpace.freeDof), ...
    'PML fine matrix is not restricted to active free DOFs.');
assert(size(pre.coarseSpace.AH, 1) == numel(pre.pml.lod.dof.coarseFree), ...
    'PML LOD coarse dimension does not match coarse free DOFs.');

relAH = norm(pre.coarseSpace.AH - pre.pml.lod.system.AH, 'fro') / ...
    max(1, norm(pre.pml.lod.system.AH, 'fro'));
assert(relAH < 1e-12, 'PML LOD coarse matrix mismatch %.3e.', relAH);

rng(7);
n = size(pre.A, 1);
r = randn(n, 2) + 1i * randn(n, 2);
y = pre.localSolver.applyInverse(r);
assert(all(isfinite(real(y(:)))) && all(isfinite(imag(y(:)))), ...
    'Local PML apply returned non-finite values.');

x = randn(n, 1) + 1i * randn(n, 1);
err = norm(pre.apply(x) - pre.applyResidual(pre.A * x)) / ...
    max(1, norm(pre.apply(x)));
assert(err < 1e-9, 'PML LXZZ identity error %.3e is too large.', err);

fprintf('PASSED: free=%d, coarse=%d, localDof=[%d,%d], relAH %.2e, identity %.2e\n', ...
    n, size(pre.coarseSpace.AH, 1), pre.local.localDofMin, ...
    pre.local.localDofMax, relAH, err);
