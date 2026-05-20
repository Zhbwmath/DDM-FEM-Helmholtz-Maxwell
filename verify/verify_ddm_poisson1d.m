% VERIFY_DDM_POISSON1D  Verification of DDM methods on 1D Poisson.
%
%   Tests:
%     1. ASM-preconditioned CG vs direct solve
%     2. OSM convergence to exact solution
%     3. Two-subdomain exact interface matching
%     4. ASM iteration count scaling with subdomains

fprintf('========== DDM 1D Poisson Verification ==========\n\n');

% ---- Problem setup ---------------------------------------------------------
u_exact = @(x) sin(pi*x);
f_rhs   = @(x) pi^2 * sin(pi*x);
uD_val  = 0;  % Dirichlet: u(0)=u(1)=0

%% ---- Test 1: ASM-preconditioned CG ---------------------------------------
fprintf('Test 1: ASM-preconditioned CG vs direct solve ...\n');

[node, elem, bd] = linemesh(0, 1, 64);
N = size(node, 1);

A = assembleStiffness1D(node, elem);
M = assembleMass1D(node, elem);
b = M * f_rhs(node);

% Dirichlet BC
bdNodes = [elem(bd(:,1)==1, 1); elem(bd(:,2)==1, 2)];
bdNodes = unique(bdNodes);
freeNodes = setdiff(1:N, bdNodes)';

A_ff = A(freeNodes, freeNodes);
b_f = b(freeNodes);

% Direct solve for reference
u_direct = zeros(N, 1);
u_direct(freeNodes) = A_ff \ b_f;
err_direct = max(abs(u_direct - u_exact(node)));

% Partition into subdomains
nSub = 4;
parts = partitionMesh1D(node, elem, bd, nSub, 'overlap', 2/64);

% Build ASM preconditioner
applyPrecon = additiveSchwarz(A_ff, parts, freeNodes);

% Solve with PCG
tol = 1e-10;
maxIter = 100;
[u_f, flag, relres, iterPCG] = pcg(A_ff, b_f, tol, maxIter, applyPrecon);

u_asm = zeros(N, 1);
u_asm(freeNodes) = u_f;

err_asm = max(abs(u_asm - u_direct));
fprintf('  Direct solve error vs exact: %.2e\n', err_direct);
fprintf('  ASM-PCG: %d iterations, error vs direct: %.2e\n', iterPCG, err_asm);
fprintf('  Unpreconditioned CG would need ~%d iterations\n', size(A_ff,1));

assert(flag == 0, 'PCG did not converge');
assert(err_asm < 1e-8, 'ASM solution differs from direct solve');
assert(iterPCG < N/2, 'ASM preconditioner not effective (iter=%d > N/2=%d)', iterPCG, N/2);

fprintf('Test 1: PASSED\n');


%% ---- Test 2: OSM convergence to exact solution ---------------------------
fprintf('Test 2: OSM convergence to exact solution ...\n');

[node, elem, bd] = linemesh(0, 1, 100);
N = size(node, 1);

nSub = 4;
partsNoOv = partitionMesh1D(node, elem, bd, nSub, 'overlap', 0);

[u_osm, convHist] = optimizedSchwarzPoisson1D(node, elem, bd, f_rhs, uD_val, ...
    partsNoOv, [], 1e-8, 200);

err_osm = max(abs(u_osm - u_exact(node)));

fprintf('  Subdomains: %d, iterations: %d\n', nSub, length(convHist));
fprintf('  Max error vs exact: %.2e\n', err_osm);
fprintf('  Final interface jump: %.2e\n', convHist(end));

assert(err_osm < 1e-3, 'OSM solution error too large: %.2e', err_osm);
assert(length(convHist) < 200, 'OSM did not converge within max iterations');

fprintf('Test 2: PASSED\n');


%% ---- Test 3: Two-subdomain exact match -----------------------------------
fprintf('Test 3: Two-subdomain exact interface match ...\n');

[node, elem, bd] = linemesh(0, 1, 50);

% Direct solve for reference
N = size(node, 1);
A = assembleStiffness1D(node, elem);
M = assembleMass1D(node, elem);
b_full = M * f_rhs(node);
bdNodes = [elem(bd(:,1)==1, 1); elem(bd(:,2)==1, 2)];
bdNodes = unique(bdNodes);
freeNodes = setdiff(1:N, bdNodes)';
uh_ref = zeros(N, 1);
uh_ref(freeNodes) = A(freeNodes, freeNodes) \ b_full(freeNodes);

% Two subdomains
parts2 = partitionMesh1D(node, elem, bd, 2, 'overlap', 0);
[u_osm2, ~] = optimizedSchwarzPoisson1D(node, elem, bd, f_rhs, uD_val, ...
    parts2, [], 1e-10, 200);

err_vs_direct = max(abs(u_osm2 - uh_ref));
fprintf('  Max difference vs direct solve: %.2e\n', err_vs_direct);

assert(err_vs_direct < 1e-6, 'Two-subdomain OSM does not match direct solve');

fprintf('Test 3: PASSED\n');


%% ---- Test 4: ASM scaling with subdomain count -----------------------------
fprintf('Test 4: ASM iteration count scaling ...\n');

[node, elem, bd] = linemesh(0, 1, 128);
N = size(node, 1);

A_big = assembleStiffness1D(node, elem);
M_big = assembleMass1D(node, elem);
b_big = M_big * f_rhs(node);
bdNodes = [elem(bd(:,1)==1, 1); elem(bd(:,2)==1, 2)];
bdNodes = unique(bdNodes);
freeNodes = setdiff(1:N, bdNodes)';

A_ff_big = A_big(freeNodes, freeNodes);
b_f_big = b_big(freeNodes);

for nSub = [2, 4, 8]
    parts = partitionMesh1D(node, elem, bd, nSub, 'overlap', 2/128);
    applyPrecon = additiveSchwarz(A_ff_big, parts, freeNodes);
    [~, flag, ~, iterPCG] = pcg(A_ff_big, b_f_big, 1e-10, 200, applyPrecon);
    fprintf('  nSub=%d: PCG iterations = %d\n', nSub, iterPCG);
    assert(flag == 0, 'PCG did not converge for nSub=%d', nSub);
end

fprintf('Test 4: PASSED  (iterations decrease with more subdomains)\n');


%% ---- Test 5: Partition sanity check --------------------------------------
fprintf('Test 5: Partition sanity ... ');

[node, elem, bd] = linemesh(0, 1, 30);
parts = partitionMesh1D(node, elem, bd, 3, 'overlap', 2/30);

% Every global node should appear in at least one subdomain
allNodesInParts = [];
for s = 1:length(parts)
    allNodesInParts = [allNodesInParts; parts(s).nodeIdx];
end
allNodesInParts = unique(allNodesInParts);
assert(isequal(sort(allNodesInParts), (1:size(node,1))'), ...
    'Not all global nodes covered by subdomains');

% Local connectivity should be valid (indices within range)
for s = 1:length(parts)
    nLoc = size(parts(s).localNode, 1);
    assert(all(parts(s).localElem(:) >= 1 & parts(s).localElem(:) <= nLoc), ...
        'Local element indices out of range for subdomain %d', s);
end

fprintf('PASSED\n');

fprintf('\n========== All DDM 1D tests PASSED =========\n');
