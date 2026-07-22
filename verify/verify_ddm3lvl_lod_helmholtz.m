% VERIFY_DDM3LVL_LOD_HELMHOLTZ  Source-faithful checks for three-level LOD-DDM Helmholtz.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

fprintf('========== Source-Faithful Three-Level LOD-DDM Helmholtz ==========\n\n');

%% ---- Shared small problem -------------------------------------------------
k = 4;
pde = helmholtzPDE(k, 'epsilon', 0, 'eta', 'sqrt');
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/4);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1/2);
lodOpts = struct('oversampling', 1, 'solveCoarse', false, ...
    'solverMode', 'direct');
partsFine = coarseHatPartition2D(node, elem, bdFlag, 1/2);
preOpts = struct('variant', 'impedance', 'coarseType', 'lod', ...
    'lodOptions', lodOpts, 'solverMode', 'direct', ...
    'adjointType', 'energy');
preExact = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, pde, ...
    partsFine, nodeH, elemH, bdH, preOpts);

sourceBase = struct('geometryMode', 'source', ...
    'sourceOversampling', lodOpts.oversampling, ...
    'tildeBufferLayers', 2, 'strictPou', true, ...
    'smax', 3, 'compareLocalBasis', true);

%% ---- Test 1: one coarse subdomain is exact -------------------------------
fprintf('Test 1: one-subdomain source geometry recovers exact coarse solve ... ');
coarseOne = buildLODCoarseSchwarzHelmholtz2D(preExact, nodeH, elemH, bdH, ...
    mergeStruct(sourceBase, struct('subdomainGrid', [1, 1])));
assert(~coarseOne.diagnostics.explicit, ...
    'Source-mode diagnostics should not explicitly assemble M0inv in this stage.');
assertSourceHierarchy(coarseOne);
rng(17);
Nc = size(coarseOne.A0, 1);
r0 = randn(Nc, 3) + 1i * randn(Nc, 3);
yExact = preExact.coarseSpace.solve(r0);
yOneSweep = coarseOne.applyM0inv(r0);
coarseOne.resetInnerSolveStats();
ySource = coarseOne.coarseSpace.solve(r0);
innerStats = coarseOne.getInnerSolveStats();
solveErr = norm(ySource - yExact, 'fro') / max(1, norm(yExact, 'fro'));
oneSweepErr = norm(yOneSweep - yExact, 'fro') / max(1, norm(yExact, 'fro'));
assert(solveErr < 1e-9, ...
    'One-subdomain source inner coarse solve is not exact: %.3e.', solveErr);
assert(oneSweepErr < 1e-9, ...
    'One-subdomain source DD preconditioner is not exact: %.3e.', oneSweepErr);
assert(strcmp(coarseOne.options.coarseSolveMode, 'oneSweep'), ...
    'Default source coarse solve mode must be oneSweep.');
assert(innerStats.calls == 0, ...
    'One-shot default should not record nested inner GMRES solves.');
fprintf('PASSED  (one-shot solve error %.3e, inner calls %d)\n', ...
    solveErr, innerStats.calls);

%% ---- Test 2: local basis recomputation is constrained --------------------
fprintf('Test 2: local LOD basis recomputation has safe-region consistency ... ');
cmp = coarseOne.basisComparison;
assert(strcmp(cmp.status, 'computed'), 'Local/global basis comparison was not computed.');
assert(cmp.maxKernelTrial < 1e-10, 'Local trial corrector violates Clement kernel.');
assert(cmp.maxKernelTest < 1e-10, 'Local test corrector violates Clement kernel.');
assert(cmp.nSafeEqualityColumns > 0, 'No safe-equality columns were classified.');
assert(cmp.maxSafeTrialRelEnergy < 1e-9, ...
    'Safe local/global trial basis mismatch is too large.');
assert(cmp.maxSafeTestRelEnergy < 1e-9, ...
    'Safe local/global test basis mismatch is too large.');
fprintf('PASSED  (safe %.3e / %.3e)\n', ...
    cmp.maxSafeTrialRelEnergy, cmp.maxSafeTestRelEnergy);

%% ---- Test 3: adjoint identity uses local source factors ------------------
fprintf('Test 3: source coarse inverse adjoint identity ... ');
x = randn(Nc, 1) + 1i * randn(Nc, 1);
y = randn(Nc, 1) + 1i * randn(Nc, 1);
adjErr = abs(x' * coarseOne.applyM0inv(y) - coarseOne.applyM0invAdjoint(x)' * y);
adjErr = adjErr / max(1, norm(x) * norm(y));
assert(adjErr < 1e-10, 'Coarse approximate inverse adjoint identity failed.');
assert(coarseOne.diagnostics.gCholFlag == 0 || ...
    coarseOne.diagnostics.gMinEigenvalue > -1e-11, ...
    'G0 is not numerically Hermitian positive semidefinite.');
fprintf('PASSED  (adjoint %.3e)\n', adjErr);

%% ---- Test 4: two-by-two source geometry and local bases ------------------
fprintf('Test 4: two-by-two source geometry satisfies cutoff hierarchy and PoU ... ');
coarseDD = buildLODCoarseSchwarzHelmholtz2D(preExact, nodeH, elemH, bdH, ...
    mergeStruct(sourceBase, struct('subdomainGrid', [2, 2], ...
    'sourceOversampling', 1, 'tildeBufferLayers', 2)));
assertSourceHierarchy(coarseDD);
assert(strcmp(coarseDD.coarseSpace.coarseSolveMode, 'oneSweep'), ...
    'Source coarseSpace.solve should use oneSweep by default.');
assert(all(strcmp({coarseDD.local.basisMode}, 'localLod')), ...
    'Source-mode local coarse matrices must use locally recomputed LOD bases.');
assert(coarseDD.basisComparison.maxKernelTrial < 1e-10, ...
    'Two-by-two local trial corrector violates Clement kernel.');
assert(coarseDD.basisComparison.maxKernelTest < 1e-10, ...
    'Two-by-two local test corrector violates Clement kernel.');
fprintf('PASSED  (safe columns %d)\n', ...
    coarseDD.basisComparison.nSafeEqualityColumns);

%% ---- Test 5: injected three-level wrapper remains algebraically coherent --
fprintf('Test 5: injected wrapper preserves LXZZ residual/function identity ... ');
threeOpts = struct('fineSpace', preExact.fineSpace, ...
    'coarseSpace', coarseDD.coarseSpace, ...
    'localSolver', preExact.localSolver, ...
    'variant', 'impedance', 'solverMode', 'direct', ...
    'adjointType', 'energy');
preThree = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, pde, ...
    partsFine, nodeH, elemH, bdH, threeOpts);
assert(isequaln(preExact.local, preThree.local), ...
    'Fine-level local solver statistics changed after coarse-solve injection.');
v = randn(size(preThree.A, 1), 1) + 1i * randn(size(preThree.A, 1), 1);
identityErr = norm(preThree.apply(v) - preThree.applyResidual(preThree.A * v)) / ...
    max(1, norm(preThree.apply(v)));
assert(identityErr < 1e-10, 'Three-level apply/applyResidual identity failed.');
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, k, 1);
coarseDD.resetInnerSolveStats();
[~, ~, relres, ~, resvec] = gmres(preThree.A, b, [], 1e-6, 20, @preThree.applyResidual);
outerInnerStats = coarseDD.getInnerSolveStats();
assert(isfinite(relres) && all(isfinite(resvec)), 'Three-level GMRES history is non-finite.');
assert(outerInnerStats.calls == 0, ...
    'One-shot source coarse solve should not run nested inner GMRES.');
fprintf('PASSED  (identity %.3e, relres %.3e, inner calls %d)\n', ...
    identityErr, relres, outerInnerStats.calls);

fprintf('\n========== Source-faithful LOD-DDM Helmholtz tests PASSED ==========\n');


function assertSourceHierarchy(method)
parts = method.parts;
Nc = size(method.A0, 1);
sumChi = zeros(Nc, 1);
for s = 1:numel(parts)
    assert(strcmp(parts(s).geometryMode, 'source'), 'Expected source geometry mode.');
    chi = parts(s).chiValues(:);
    chiGreater = parts(s).chiGreaterValues(:);
    sumChi = sumChi + chi;
    suppChi = chi > 0;
    plateauGreater = chiGreater == 1;
    assert(all(suppChi <= plateauGreater), ...
        'supp chi is not contained in {chiGreater = 1}.');
    assert(all(plateauGreater <= parts(s).omegaNodes(:)), ...
        '{chiGreater = 1} is not contained in Omega_0,l.');
    assert(all(parts(s).omegaNodes(:) <= parts(s).tildeNodes(:)), ...
        'Omega_0,l is not contained in tildeOmega_0,l.');
    assert(parts(s).tildeLayers == ...
        method.options.sourceOversampling + method.options.tildeBufferLayers, ...
        'tildeOmega_0,l was not expanded by m + C coarse layers.');
end
assert(all(sumChi == 1), ...
    'sum_l Pi_H(chi_l .) is not exactly the identity on coarse coefficients.');
assert(method.stats.pouExact, 'Source-mode PoU exactness flag was not set.');
end


function out = mergeStruct(a, b)
out = a;
names = fieldnames(b);
for i = 1:numel(names)
    out.(names{i}) = b.(names{i});
end
end
