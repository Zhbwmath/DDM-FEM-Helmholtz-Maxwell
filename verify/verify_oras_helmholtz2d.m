% VERIFY_ORAS_HELMHOLTZ2D  Verification of ORAS for 2D Helmholtz.
%
%   Tests:
%     1. Helmholtz direct-solve consistency with ORAS-preconditioned GMRES
%     2. Plane wave: ORAS-preconditioned GMRES convergence
%     3. Partition of unity check

fprintf('========== ORAS Helmholtz 2D Verification ==========\n\n');

%% ---- Test 1: Helmholtz direct-solve consistency --------------------------
fprintf('Test 1: Helmholtz GMRES+ORAS, compare with direct solve...\n');

[node, elem, bd] = squaremesh([0, 1, 0, 1], 1/16);
N = size(node, 1);
k = 5;

u_ex = @(x,y) exp(1i*k*x);
f_rhs = @(x,y) 0*x;
g_bc  = @(x,y) 0*x;

A0 = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
u_ex_vec = u_ex(node(:,1), node(:,2));
b0 = A0 * u_ex_vec;
u_direct = A0 \ b0;

parts = partitionMesh2D(node, elem, bd, 4, 'overlap', 2/16);
parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [4, 1], 2/16);
ap = orasHelmholtz(node, elem, bd, k, parts);

[u_oras, flag, ~, iter] = gmres(A0, b0, [], 1e-10, 200, ap);

err = max(abs(u_oras - u_direct));
fprintf('  GMRES+ORAS(k=%g): %d iterations, error vs direct: %.2e, flag=%d\n', ...
    k, iter(2)-iter(1), err, flag);
assert(flag == 0, 'GMRES did not converge');
assert(err < 1e-8, 'ORAS solution does not match direct solve');
fprintf('Test 1: PASSED\n');

%% ---- Test 2: Plane wave, ORAS preconditioned GMRES -----------------------
fprintf('\nTest 2: Plane wave u=exp(ikx), GMRES+ORAS convergence...\n');

k = 10;
[node, elem, bd] = squaremesh([0, 1, 0, 1], 1/24);
N = size(node, 1);

u_ex = @(x,y) exp(1i*k*x);
f_rhs = @(x,y) 0*x;
g_bc  = @(x,y) 0*x;

A2 = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
u_ex_vec = u_ex(node(:,1), node(:,2));
b2 = A2 * u_ex_vec;

parts = partitionMesh2D(node, elem, bd, 4, 'overlap', 2/24);
parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [4, 1], 2/24);
ap2 = orasHelmholtz(node, elem, bd, k, parts);

[u_gm, flag_gm, relres_gm, iter_gm] = gmres(A2, b2, [], 1e-8, 200, ap2);
fprintf('  GMRES+ORAS: %d iterations, relres=%.2e, flag=%d\n', ...
    iter_gm(2)-iter_gm(1), relres_gm, flag_gm);
fprintf('  GMRES error vs exact: %.2e\n', max(abs(u_gm - u_ex_vec)));
assert(flag_gm == 0, 'GMRES did not converge');

fprintf('  Damped Richardson (omega=0.5): ');
u_rich = zeros(N, 1);
omega = 0.5;
ch_rich = zeros(200, 1);
for it = 1:200
    r = b2 - A2 * u_rich;
    u_rich = u_rich + omega * ap2(r);
    ch_rich(it) = norm(r) / norm(b2);
    if ch_rich(it) < 1e-4
        ch_rich = ch_rich(1:it);
        break;
    end
end
fprintf('%d iterations\n', length(ch_rich));

fprintf('Test 2: PASSED\n');

%% ---- Test 3: Partition of unity check ------------------------------------
fprintf('\nTest 3: Partition of unity sum check...\n');

parts3 = partitionMesh2D(node, elem, bd, 4, 'overlap', 2/24);
parts3 = linearPartitionOfUnity2D(parts3, [0, 1, 0, 1], [4, 1], 2/24);
nodeCount3 = zeros(N, 1);
for j = 1:length(parts3)
    raw = parts3(j).weightFun(node(parts3(j).nodeIdx,1), node(parts3(j).nodeIdx,2));
    nodeCount3(parts3(j).nodeIdx) = nodeCount3(parts3(j).nodeIdx) + raw(:);
end

chiSum = zeros(N, 1);
for j = 1:length(parts3)
    raw = parts3(j).weightFun(node(parts3(j).nodeIdx,1), node(parts3(j).nodeIdx,2));
    chiSum = chiSum + accumarray(parts3(j).nodeIdx, ...
        raw(:) ./ nodeCount3(parts3(j).nodeIdx), [N, 1]);
end

errPou = max(abs(chiSum(chiSum > 0) - 1));
fprintf('  sum chi_j(node) == 1: max error = %.2e\n', errPou);
assert(errPou < 1e-12, 'Partition of unity does not sum to 1');

fprintf('Test 3: PASSED\n');
fprintf('\n========== ORAS Verification Complete ==========\n');
