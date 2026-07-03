% VERIFY_DDM3LVL_LOD_HELMHOLTZ  Smoke checks for three-level LOD-DDM Helmholtz.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

fprintf('========== Three-Level LOD-DDM Helmholtz Smoke ==========\n\n');

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

%% ---- Test 1: one coarse subdomain is exact -------------------------------
fprintf('Test 1: one-subdomain coarse Schwarz recovers exact coarse solve ... ');
coarseOpts = struct('subdomainGrid', [1, 1], 'smax', 3, ...
    'compareLocalBasis', true);
coarseOne = buildLODCoarseSchwarzHelmholtz2D(preExact, nodeH, elemH, bdH, coarseOpts);
assert(coarseOne.diagnostics.explicit, 'Expected explicit diagnostics on smoke case.');
eNorm = norm(full(coarseOne.diagnostics.E0), 'fro');
assert(eNorm < 1e-9, 'One-subdomain E0 is not numerically zero: %.3e.', eNorm);
assert(coarseOne.diagnostics.sContract == 1, 'One-subdomain contraction should occur at s=1.');
fprintf('PASSED  (||E0||_F %.3e)\n', eNorm);

%% ---- Test 2: local basis recomputation is constrained --------------------
fprintf('Test 2: local LOD basis recomputation satisfies local kernel constraints ... ');
cmp = coarseOne.basisComparison;
assert(strcmp(cmp.status, 'computed'), 'Local/global basis comparison was not computed.');
assert(cmp.maxKernelTrial < 1e-10, 'Local trial corrector violates Clement kernel.');
assert(cmp.maxKernelTest < 1e-10, 'Local test corrector violates Clement kernel.');
assert(cmp.maxContainedTrialRelEnergy < 1e-9, ...
    'Contained local/global trial basis mismatch is too large.');
assert(cmp.maxContainedTestRelEnergy < 1e-9, ...
    'Contained local/global test basis mismatch is too large.');
fprintf('PASSED  (kernel %.3e / %.3e)\n', cmp.maxKernelTrial, cmp.maxKernelTest);

%% ---- Test 3: adjoint and s-sweep identities ------------------------------
fprintf('Test 3: coarse inverse adjoint and s-sweep identities ... ');
rng(31);
Nc = size(coarseOne.A0, 1);
x = randn(Nc, 1) + 1i * randn(Nc, 1);
y = randn(Nc, 1) + 1i * randn(Nc, 1);
adjErr = abs(x' * coarseOne.applyM0inv(y) - coarseOne.applyM0invAdjoint(x)' * y);
adjErr = adjErr / max(1, norm(x) * norm(y));
G2 = explicitOperator(@(r) coarseOne.applyG0s(r, 2), Nc);
idErr = norm(speye(Nc) - G2 * coarseOne.A0 - coarseOne.diagnostics.E0^2, 'fro');
assert(adjErr < 1e-10, 'Coarse approximate inverse adjoint identity failed.');
assert(idErr < 1e-10, 's-sweep residual identity failed.');
assert(coarseOne.diagnostics.gMinEigenvalue > -1e-11, 'G0 is not numerically SPD.');
fprintf('PASSED  (adjoint %.3e, identity %.3e)\n', adjErr, idErr);

%% ---- Test 4: nontrivial coarse DD diagnostics are finite -----------------
fprintf('Test 4: two-by-two coarse Schwarz diagnostics are finite ... ');
coarseOpts = struct('subdomainGrid', [2, 2], 'overlap', 0.25, ...
    'greaterOverlap', 0.25, 'smax', 4, 'compareLocalBasis', true);
coarseDD = buildLODCoarseSchwarzHelmholtz2D(preExact, nodeH, elemH, bdH, coarseOpts);
assert(all(isfinite(coarseDD.diagnostics.normEPower)), 'Power norms contain non-finite values.');
assert(all(isfinite(coarseDD.diagnostics.alpha)), 'FOV alpha values contain non-finite values.');
assert(coarseDD.basisComparison.maxKernelTrial < 1e-10, ...
    'Two-by-two local trial corrector violates Clement kernel.');
assert(coarseDD.basisComparison.maxKernelTest < 1e-10, ...
    'Two-by-two local test corrector violates Clement kernel.');
fprintf('PASSED  (s0=%s, alpha>0 at %s)\n', ...
    numberOrDash(coarseDD.diagnostics.sContract), ...
    numberOrDash(coarseDD.diagnostics.sFovPositive));

%% ---- Test 5: injected three-level wrapper remains algebraically coherent --
fprintf('Test 5: injected three-level wrapper preserves LXZZ residual/function identity ... ');
threeOpts = struct('fineSpace', preExact.fineSpace, ...
    'coarseSpace', coarseDD.coarseSpace, ...
    'localSolver', preExact.localSolver, ...
    'variant', 'impedance', 'solverMode', 'direct', ...
    'adjointType', 'energy');
preThree = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, pde, ...
    partsFine, nodeH, elemH, bdH, threeOpts);
v = randn(size(preThree.A, 1), 1) + 1i * randn(size(preThree.A, 1), 1);
identityErr = norm(preThree.apply(v) - preThree.applyResidual(preThree.A * v)) / ...
    max(1, norm(preThree.apply(v)));
assert(identityErr < 1e-10, 'Three-level apply/applyResidual identity failed.');
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, k, 1);
[~, ~, relres, ~, resvec] = gmres(preThree.A, b, [], 1e-6, 20, @preThree.applyResidual);
assert(isfinite(relres) && all(isfinite(resvec)), 'Three-level GMRES history is non-finite.');
fprintf('PASSED  (identity %.3e, relres %.3e)\n', identityErr, relres);

fprintf('\n========== Three-Level LOD-DDM Helmholtz smoke tests PASSED ==========\n');


function B = explicitOperator(applyFun, n)
B = applyFun(speye(n));
B = sparse(B);
end


function s = numberOrDash(x)
if isnan(x)
    s = '-';
else
    s = sprintf('%d', x);
end
end
