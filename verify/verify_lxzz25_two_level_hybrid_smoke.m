% VERIFY_LXZZ25_TWO_LEVEL_HYBRID_SMOKE  Smoke checks for LXZZ25 hybrid LOD-DDM.

fprintf('========== LXZZ25 Two-Level Hybrid Schwarz Smoke ==========\n\n');

%% ---- Test 1: Matrix assembly --------------------------------------------
fprintf('Test 1: Helmholtz assembly matrix properties ... ');
k = 4;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/8);
K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
Mb = assembleBoundaryMass2D(node, elem, bdFlag, 1);
A = K - k^2 * M - 1i * k * Mb;

assert(norm(K - K', 'fro') < 1e-12 * max(1, norm(K, 'fro')), 'K is not Hermitian.');
assert(norm(M - M', 'fro') < 1e-12 * max(1, norm(M, 'fro')), 'M is not Hermitian.');
assert(norm(Mb - Mb', 'fro') < 1e-12 * max(1, norm(Mb, 'fro')), 'Mb is not Hermitian.');
assert(min(eig(full(M))) > 0, 'M is not positive definite on this mesh.');
assert(min(eig(full(K + 1e-12 * speye(size(K))))) > -1e-10, 'K is not semidefinite.');
assert(min(eig(full(Mb + 1e-12 * speye(size(Mb))))) > -1e-10, 'Mb is not semidefinite.');
assert(norm(A - A', 'fro') > 1e-8, 'Impedance Helmholtz matrix should be non-Hermitian.');
fprintf('PASSED\n');

%% ---- Test 2: Plane-wave consistency --------------------------------------
fprintf('Test 2: Plane-wave consistency under refinement ... ');
err = zeros(2, 1);
hs = [1/8, 1/16];
for i = 1:numel(hs)
    [nodei, elemi, bdi] = squaremesh([0, 1, 0, 1], hs(i));
    Ai = assembleHelmholtz2D(nodei, elemi, bdi, k, 0, [], 1);
    bi = assemblePlaneWaveBoundaryLoadP1(nodei, elemi, bdi, k);
    uh = Ai \ bi;
    ue = planeWave(nodei, k);
    Ei = assembleStiffness2D(nodei, elemi, 1) + k^2 * assembleMass2D(nodei, elemi, 1);
    err(i) = sqrt(real((uh - ue)' * Ei * (uh - ue))) / sqrt(real(ue' * Ei * ue));
end
assert(all(isfinite(err)), 'Plane-wave errors are not finite.');
assert(err(2) < 0.9 * err(1), 'Plane-wave error did not decrease enough under refinement.');
fprintf('PASSED  (%.3e -> %.3e)\n', err(1), err(2));

%% ---- Test 3: LOD interpolation kernel ------------------------------------
fprintf('Test 3: LOD correctors satisfy Clement kernel constraints ... ');
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1/2);
lodOpts = struct('oversampling', 1, 'solveCoarse', false, 'solverMode', 'direct');
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, node, elem, bdFlag, k, 0, 0, lodOpts);
P = prolongateNestedP1(nodeH, elemH, node);
Q = weightedClementP1(node, elem, nodeH, elemH);
trialCorrector = P - lod.basis.trial;
testCorrector = P - lod.basis.test;
kernelResidual = max(norm(Q * trialCorrector, 'fro'), norm(Q * testCorrector, 'fro'));
assert(kernelResidual < 1e-10, 'LOD corrector is not in the Clement kernel.');
fprintf('PASSED  (kernel residual %.3e)\n', kernelResidual);

%% ---- Test 4: Patch monotonicity ------------------------------------------
fprintf('Test 4: LOD patches enlarge monotonically ... ');
p0 = lodBuildPatches2D(nodeH, elemH, node, elem, bdFlag, 0);
p1 = lodBuildPatches2D(nodeH, elemH, node, elem, bdFlag, 1);
p2 = lodBuildPatches2D(nodeH, elemH, node, elem, bdFlag, 2);
for T = 1:size(elemH, 1)
    assert(numel(p0.fineElemIds{T}) <= numel(p1.fineElemIds{T}), 'ell=1 patch shrank.');
    assert(numel(p1.fineElemIds{T}) <= numel(p2.fineElemIds{T}), 'ell=2 patch shrank.');
end
fprintf('PASSED\n');

%% ---- Test 5: Hybrid wrapper sanity ---------------------------------------
fprintf('Test 5: Corrected residual and function-level hybrid identities ... ');
parts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', 1/8);
parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [2, 2], 1/8);

opts = struct('lod', lod, 'lodOptions', lodOpts, 'solverMode', 'lu');
opts.variant = 'dirichlet';
q1 = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, parts, nodeH, elemH, bdH, opts);
opts.variant = 'impedance';
q2 = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, parts, nodeH, elemH, bdH, opts);

rng(1);
x = randn(size(node, 1), 1) + 1i * randn(size(node, 1), 1);
y = randn(size(node, 1), 1) + 1i * randn(size(node, 1), 1);
for q = {q1, q2}
    z = q{1}.apply(x);
    zr = q{1}.applyResidual(q{1}.A * x);
    linErr = norm(q{1}.apply(x + y) - q{1}.apply(x) - q{1}.apply(y)) / max(1, norm(z));
    adjErr = abs(((x - q{1}.applyQ0(x))' * q{1}.energy * y) ...
        - (x' * q{1}.energy * q{1}.applyEnergyAdjointIMinusQ0(y)));
    adjErr = adjErr / max(1, norm(x) * norm(y));
    z0 = q{1}.applyM0Inverse(y);
    rc = y - q{1}.A * z0;
    wrongRc = y - q{1}.applyQ0(y);
    assert(all(isfinite(z)), 'Hybrid wrapper produced non-finite entries.');
    assert(norm(z - zr) < 1e-10 * max(1, norm(z)), 'B_i^{-1} A is inconsistent with Q_m^{(i)}.');
    assert(norm(rc - wrongRc) > 1e-8 * max(1, norm(rc)), 'Residual update did not distinguish M0^{-1} from Q0.');
    assert(linErr < 1e-10, 'Hybrid wrapper is not linear.');
    assert(adjErr < 1e-8, '(I-Q0) D_kappa-adjoint identity failed.');
end
fprintf('PASSED\n');

%% ---- Test 6: Small explicit preconditioned operators ---------------------
fprintf('Test 6: Explicit small matrices satisfy B_i^{-1} A = Q_m^{(i)} ... ');
for q = {q1, q2}
    BInv = explicitOperator(q{1}.applyResidual, size(node, 1));
    Qm = explicitOperator(q{1}.apply, size(node, 1));
    identityErr = norm(BInv * q{1}.A - Qm, 'fro') / max(1, norm(Qm, 'fro'));
    assert(all(isfinite(nonzeros(BInv))), 'Explicit residual inverse matrix has non-finite entries.');
    assert(all(isfinite(nonzeros(Qm))), 'Explicit hybrid operator matrix has non-finite entries.');
    assert(identityErr < 1e-10, 'Explicit identity B_i^{-1} A = Q_m^{(i)} failed.');
end
fprintf('PASSED\n');

fprintf('\n========== LXZZ25 hybrid smoke tests PASSED ==========\n');


function u = planeWave(node, k)
d = [1 / sqrt(2), 1 / sqrt(2)];
u = exp(1i * k * (node(:, 1) * d(1) + node(:, 2) * d(2)));
end


function b = assemblePlaneWaveBoundaryLoadP1(node, elem, bdFlag, k)
N = size(node, 1);
edgeVertex = [2 3; 3 1; 1 2];
normals = [1, 0; 0, -1; -1, 0; 0, 1];
b = zeros(N, 1);
for t = 1:size(elem, 1)
    for e = 1:3
        if bdFlag(t, e) ~= 1, continue; end
        va = elem(t, edgeVertex(e, 1));
        vb = elem(t, edgeVertex(e, 2));
        mid = 0.5 * (node(va, :) + node(vb, :));
        normal = squareBoundaryNormal(mid, normals);
        ga = planeWaveBoundaryValue(node(va, :), normal, k);
        gb = planeWaveBoundaryValue(node(vb, :), normal, k);
        L = norm(node(vb, :) - node(va, :));
        b([va; vb]) = b([va; vb]) + L / 6 * [2, 1; 1, 2] * [ga; gb];
    end
end
end


function n = squareBoundaryNormal(x, normals)
tol = 1e-12;
if abs(x(1) - 1) < tol
    n = normals(1, :);
elseif abs(x(2)) < tol
    n = normals(2, :);
elseif abs(x(1)) < tol
    n = normals(3, :);
elseif abs(x(2) - 1) < tol
    n = normals(4, :);
else
    error('verify_lxzz25:normal', 'Boundary midpoint was not on the unit square.');
end
end


function g = planeWaveBoundaryValue(x, n, k)
d = [1 / sqrt(2), 1 / sqrt(2)];
u = exp(1i * k * (x(:, 1) * d(1) + x(:, 2) * d(2)));
g = 1i * k * (d * n.' - 1) .* u;
end


function B = explicitOperator(applyFun, n)
B = zeros(n, n);
for j = 1:n
    ej = zeros(n, 1);
    ej(j) = 1;
    B(:, j) = applyFun(ej);
end
B = sparse(B);
end
