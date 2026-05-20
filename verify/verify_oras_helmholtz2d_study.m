% VERIFY_ORAS_HELMHOLTZ2D_STUDY  ORAS/GMRES comparison for Helmholtz.

fprintf('========== ORAS Helmholtz 2D Study ==========\n\n');
fprintf('References: ORAS uses local impedance conditions; compare overlap, k, and partition shape.\n');

f_rhs = @(x, y) 0*x;
g_bc = @(x, y) 0*x;
gmresTol = 1e-8;
maxIter = 200;

%% Direct consistency
fprintf('Test 1: direct consistency, k=10, h=1/24, checkerboard 2x2, delta=2h ... ');
k = 10;
h = 1/24;
[node, elem, bd] = squaremesh([0, 1, 0, 1], h);
A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
uEx = exp(1i*k*node(:,1));
b = A * uEx;
uDirect = A \ b;
parts = partitionMesh2D(node, elem, bd, [2, 2], 'overlap', 2*h);
parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [2, 2], 2*h);
orasPrec = orasHelmholtz(node, elem, bd, k, parts);
[uGM, flag, relres, iter] = gmres(A, b, [], gmresTol, maxIter, orasPrec);
assert(flag == 0, 'ORAS-GMRES failed direct consistency, relres %.3e', relres);
assert(norm(uGM-uDirect)/norm(uDirect) < 1e-8, 'ORAS-GMRES differs from direct solve.');
fprintf('PASSED (%d iterations)\n', iter(2));

%% Overlap comparison
fprintf('\nTable 1: ORAS overlap comparison, k=20, h=1/32, checkerboard 2x2\n');
k = 20;
h = 1/32;
[node, elem, bd] = squaremesh([0, 1, 0, 1], h);
A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
uEx = exp(1i*k*node(:,1));
b = A * uEx;
fprintf('%-8s %-10s %-10s %-10s\n', 'delta/h', 'GMRES_it', 'relres', 'flag');
fprintf('%s\n', repmat('-', 1, 45));
for m = [0, 1, 2, 4]
    parts = partitionMesh2D(node, elem, bd, [2, 2], 'overlap', m*h);
    parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [2, 2], m*h);
    orasPrec = orasHelmholtz(node, elem, bd, k, parts);
    [~, flag, relres, iter] = gmres(A, b, [], gmresTol, maxIter, orasPrec);
    fprintf('%-8d %-10d %-10.2e %-10d\n', m, iter(2), relres, flag);
    assert(flag == 0, 'ORAS overlap case delta/h=%d failed.', m);
end

%% k comparison at fixed points per wavelength
fprintf('\nTable 2: ORAS k comparison, about 10 points per wavelength, checkerboard 2x2, delta=2h\n');
fprintf('%-8s %-8s %-10s %-10s\n', 'k', '1/h', 'N', 'GMRES_it');
fprintf('%s\n', repmat('-', 1, 42));
for k = [10, 20, 30]
    h = min(1/24, 1/(10*k/(2*pi)));
    nSide = ceil(1/h);
    h = 1/nSide;
    [node, elem, bd] = squaremesh([0, 1, 0, 1], h);
    A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
    uEx = exp(1i*k*node(:,1));
    b = A * uEx;
    parts = partitionMesh2D(node, elem, bd, [2, 2], 'overlap', 2*h);
    parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [2, 2], 2*h);
    orasPrec = orasHelmholtz(node, elem, bd, k, parts);
    [~, flag, relres, iter] = gmres(A, b, [], gmresTol, maxIter, orasPrec);
    fprintf('%-8d %-8d %-10d %-10d\n', k, nSide, size(node, 1), iter(2));
    assert(flag == 0, 'ORAS k=%d failed, relres %.3e', k, relres);
end

%% Partition shape comparison
fprintf('\nTable 3: ORAS partition comparison, k=20, h=1/32, delta=2h\n');
k = 20;
h = 1/32;
[node, elem, bd] = squaremesh([0, 1, 0, 1], h);
A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
uEx = exp(1i*k*node(:,1));
b = A * uEx;
configs = {
    'strip-4', 4
    'grid-2x2', [2, 2]
    'grid-4x2', [4, 2]
};
fprintf('%-12s %-10s %-10s\n', 'partition', 'subdomains', 'GMRES_it');
fprintf('%s\n', repmat('-', 1, 40));
for i = 1:size(configs, 1)
    parts = partitionMesh2D(node, elem, bd, configs{i, 2}, 'overlap', 2*h);
    if isscalar(configs{i, 2})
        gridSize = [configs{i, 2}, 1];
    else
        gridSize = configs{i, 2};
    end
    parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], gridSize, 2*h);
    orasPrec = orasHelmholtz(node, elem, bd, k, parts);
    [~, flag, relres, iter] = gmres(A, b, [], gmresTol, maxIter, orasPrec);
    fprintf('%-12s %-10d %-10d\n', configs{i, 1}, length(parts), iter(2));
    assert(flag == 0, 'ORAS partition %s failed, relres %.3e', configs{i, 1}, relres);
end

fprintf('\n========== ORAS Helmholtz study PASSED ==========\n');
