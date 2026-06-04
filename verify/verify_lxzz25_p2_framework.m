% VERIFY_LXZZ25_P2_FRAMEWORK  Check abstract LXZZ wrapper on P2 fine space.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

fprintf('========== LXZZ25 P2 fine-space framework verification ==========\n\n');

k = envNumber('LXZZ25_P2_SMOKE_K', 8);
hInv = alignFineInv(ceil(k^(3/2)), k);
HInv = k;
m = envNumber('LXZZ25_P2_SMOKE_M', 2);

[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1 / hInv);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1 / HInv);
[node2, elem2] = extendMesh2D(node, elem, 2);

fprintf('Base mesh: k=%g h=1/%d H=1/%d P2 dofs=%d\n', ...
    k, hInv, HInv, size(node2, 1));

%% ---- Test 1: P2 matrices and P1-to-P2 embedding -------------------------
fprintf('Test 1: P2 operator and P1-to-P2 embedding ... ');
K2 = assembleStiffness2D(node2, elem2, 2);
M2 = assembleMass2D(node2, elem2, 2);
Mb2 = assembleBoundaryMass2D(node2, elem2, bdFlag, 2);
A2 = K2 - k^2 * M2 - 1i * k * Mb2;
E21 = prolongate_P1_P2(node, elem);
u1 = 1 + node(:, 1) - 2 * node(:, 2);
u2 = E21 * u1;
u2Exact = 1 + node2(:, 1) - 2 * node2(:, 2);
assert(size(A2, 1) == size(node2, 1), 'P2 Helmholtz matrix has wrong size.');
assert(all(isfinite(nonzeros(A2))), 'P2 Helmholtz matrix contains non-finite entries.');
assert(norm(u2 - u2Exact, inf) < 1e-13, 'P1-to-P2 embedding is not exact for linear functions.');
fprintf('PASSED\n');

%% ---- Test 2: Embedded LOD coarse matrix is recomputed in P2 --------------
fprintf('Test 2: embedded LOD basis and P2 coarse matrix ... ');
parts = coarseHatPartition2D(node, elem, bdFlag, 1 / HInv);
lodOpts = struct('oversampling', m, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', false);
opts = struct('degree', 2, 'variant', 'dirichlet', 'coarseType', 'lod', ...
    'lodOptions', lodOpts, 'solverMode', 'adaptive', 'useParfor', false, ...
    'adjointType', 'reference');
preD = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, opts);
AHcheck = preD.basis.test' * preD.A * preD.basis.trial;
relAH = norm(AHcheck - preD.basis.AH, 'fro') / max(1, norm(preD.basis.AH, 'fro'));
assert(size(preD.basis.trial, 1) == size(node2, 1), 'Embedded LOD trial basis is not in P2 space.');
assert(relAH < 1e-12, 'P2 coarse matrix was not recomputed from embedded basis.');
fprintf('PASSED  (relAH %.3e)\n', relAH);

%% ---- Test 3: P2 Dirichlet masks and impedance partition weights ----------
fprintf('Test 3: P2 local DOF masks and weights ... ');
dirCounts = zeros(numel(parts), 1);
weightSum = zeros(size(node2, 1), 1);
covered = false(size(node2, 1), 1);
for s = 1:numel(parts)
    idx = unique(elem2(parts(s).elemIdx, :));
    w = max(parts(s).weightFun(node2(idx, 1), node2(idx, 2)), 0);
    dirCounts(s) = nnz(w > 1e-12);
    assert(all(w(w > 1e-12) > 1e-12), 'Dirichlet positive-hat mask kept zero-weight DOFs.');
    weightSum(idx) = weightSum(idx) + w(:);
    covered(idx) = true;
end
assert(abs(min(dirCounts) - preD.local.localDofMin) <= 0, ...
    'Dirichlet local DOF count does not match positive-hat mask.');
assert(all(weightSum(covered) > 0), 'Some covered P2 DOFs have zero accumulated impedance weight.');
fprintf('PASSED\n');

%% ---- Test 4: P2 Q1/Q2 small GMRES smoke ---------------------------------
fprintf('Test 4: P2 Q1/Q2 Euclidean GMRES smoke ...\n');
b2 = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, k, 2);
variants = {'dirichlet', 'impedance'};
for i = 1:numel(variants)
    opts.variant = variants{i};
    pre = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, ...
        parts, nodeH, elemH, bdH, opts);
    tApply = tic;
    [~, flag, relres, iter] = gmres(pre.A, b2, [], 1e-6, 30, @pre.applyResidual);
    applyS = toc(tApply);
    iterCount = gmresIterationCount(iter);
    fprintf('  %-10s flag=%d iter=%d relres=%.3e setup=(%.2fs, %.2fs) solve=%.2fs localMode=%s\n', ...
        variants{i}, flag, iterCount, relres, pre.timing.coarseSetup, ...
        pre.timing.localSetup, applyS, pre.local.solverModeEffective);
    assert(flag == 0 || relres < 1e-5, 'P2 %s GMRES smoke did not converge enough.', variants{i});
end

fprintf('\n========== LXZZ25 P2 fine-space framework tests PASSED ==========\n');


function hInv = alignFineInv(raw, divisors)
hInv = max(1, ceil(raw));
for d = divisors(:).'
    if d <= 0, continue; end
    hInv = ceil(hInv / d) * d;
end
end


function n = envNumber(name, defaultValue)
txt = getenv(name);
if isempty(txt)
    n = defaultValue;
else
    n = str2double(txt);
    if isnan(n), n = defaultValue; end
end
end


function n = gmresIterationCount(iter)
if numel(iter) == 2
    n = iter(2);
else
    n = iter;
end
end
