% VERIFY_ORAS_STUDY  Comprehensive ORAS parameter study for 2D Helmholtz.
%
%   -(Δ + k²)u = f,  ∂u/∂n - iku = g on ∂Ω
%   ORAS preconditioner + GMRES (as in Gong-Graham-Spence 2022)
%
%   Studies:
%     1. h-refinement: GMRES iters vs h (fixed k, H, δ)
%     2. k-dependence: GMRES iters vs k (fixed h, H, δ)
%     3. Overlap δ effect: GMRES iters vs δ/h
%     4. Subdomain size H effect
%     5. Strip vs Checkerboard comparison

fprintf('========== ORAS Helmholtz 2D Parameter Study ==========\n');
fprintf('Equation: -(Delta + k^2)u = 0,  d/dn - iku = g\n');
fprintf('Exact: u = exp(ikx), plane wave in +x direction\n\n');

u_ex = @(x,y) exp(1i*0*x);  % placeholder
f_rhs = @(x,y) 0*x;
g_bc  = @(x,y) 0*x;
gmresTol = 1e-8;  maxIter = 300;

%% Table 1: h-refinement — GMRES iters vs h (fixed k, H, δ) ----------------
fprintf('==============================================================\n');
fprintf('TABLE 1: GMRES iterations vs h (k=20, strip 4, delta=2h)\n');
fprintf('==============================================================\n');

k = 20;
fprintf('%-6s %-6s %-6s %-5s %-6s %-6s\n', '1/h','N','nSub','H','delta','GMRES');
fprintf('%s\n', repmat('-', 1, 45));

for nSide = [16, 24, 32, 48]
    h = 1 / nSide;
    delta = 2 * h;  % δ = 2h per rule
    [node, elem, bd] = squaremesh([0,1,0,1], h);
    N = size(node, 1);

    A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
    u_ex_vec = u_ex(node(:,1), node(:,2));
    u_ex_vec = exp(1i*k*node(:,1));  % u = exp(ikx)
    b = A * u_ex_vec;

    parts = partitionMesh2D(node, elem, bd, 4, 'overlap', delta);
    ap = orasHelmholtz(node, elem, bd, k, parts);

    [~, flag, ~, iter] = gmres(A, b, [], gmresTol, maxIter, ap);
    its = iter(2) - iter(1);
    fprintf('%-6d %-6d %-5d %-5.2f %-6.4f %-6d\n', nSide, N, 4, 1/4, delta, its);
    assert(flag == 0, 'GMRES did not converge');
end

%% Table 2: k-dependence — GMRES iters vs k (fixed h, H, δ) -----------------
fprintf('\n==============================================================\n');
fprintf('TABLE 2: GMRES iterations vs k (h=1/32, strip 4, delta=2h)\n');
fprintf('==============================================================\n');

nSide = 32;  h = 1/nSide;  delta = 2*h;
[node, elem, bd] = squaremesh([0,1,0,1], h);
N = size(node, 1);
parts = partitionMesh2D(node, elem, bd, 4, 'overlap', delta);

fprintf('%-6s %-6s %-6s %-6s %-6s\n', 'k','N','nSub','delta','GMRES');
fprintf('%s\n', repmat('-', 1, 40));

for k = [5, 10, 20, 30, 40, 50]
    A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
    u_ex_vec = exp(1i*k*node(:,1));
    b = A * u_ex_vec;
    ap = orasHelmholtz(node, elem, bd, k, parts);
    [~, flag, ~, iter] = gmres(A, b, [], gmresTol, maxIter, ap);
    its = iter(2) - iter(1);
    fprintf('%-6d %-6d %-5d %-6.4f %-6d\n', k, N, 4, delta, its);
    assert(flag == 0, 'GMRES did not converge');
end

%% Table 3: Overlap delta effect — GMRES vs δ/h ----------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 3: GMRES iterations vs overlap delta (k=20, h=1/32, strip 4)\n');
fprintf('==============================================================\n');

k = 20;
fprintf('%-8s %-6s %-6s %-6s\n', 'delta/h','delta','GMRES','cond_est');
fprintf('%s\n', repmat('-', 1, 40));

for m = [0, 1, 2, 3, 4, 6]
    delta = m * h;
    if delta == 0, delta = 0; end  % exact zero
    parts_d = partitionMesh2D(node, elem, bd, 4, 'overlap', delta);
    A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
    u_ex_vec = exp(1i*k*node(:,1));
    b = A * u_ex_vec;
    ap_d = orasHelmholtz(node, elem, bd, k, parts_d);
    [~, flag, ~, iter] = gmres(A, b, [], gmresTol, maxIter, ap_d);
    its = iter(2) - iter(1);
    % Estimate condition number from GMRES convergence
    condEst = condest(A);  % rough estimate
    fprintf('%-8d %-6.4f %-6d %-6.0f\n', m, delta, its, condEst);
end

%% Table 4: Subdomain size H effect -----------------------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 4: GMRES iterations vs subdomain count (k=20, h=1/32, delta=2h)\n');
fprintf('==============================================================\n');

delta = 2*h;
fprintf('%-5s %-6s %-5s %-6s\n', 'nSub','H','delta','GMRES');
fprintf('%s\n', repmat('-', 1, 35));

for nSub = [2, 3, 4, 6, 8]
    H = 1 / nSub;
    parts_h = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);
    A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
    u_ex_vec = exp(1i*k*node(:,1));
    b = A * u_ex_vec;
    ap_h = orasHelmholtz(node, elem, bd, k, parts_h);
    [~, flag, ~, iter] = gmres(A, b, [], gmresTol, maxIter, ap_h);
    its = iter(2) - iter(1);
    fprintf('%-5d %-6.2f %-5.4f %-6d\n', nSub, H, delta, its);
end

%% Table 5: Strip vs Checkerboard -------------------------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 5: Strip vs Checkerboard (k=20, h=1/32, delta=2h)\n');
fprintf('==============================================================\n');

delta = 2*h;
fprintf('%-14s %-5s %-5s %-6s\n', 'Partition','nSub','H_eff','GMRES');
fprintf('%s\n', repmat('-', 1, 40));

for cfg = { {'Strip 4', 4}, {'Grid 2x2', [2,2]}, {'Grid 3x3', [3,3]} }
    name = cfg{1}{1};  nSub = cfg{1}{2};
    parts_c = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);
    A = assembleHelmholtz2D(node, elem, bd, k, f_rhs, g_bc);
    u_ex_vec = exp(1i*k*node(:,1));
    b = A * u_ex_vec;
    ap_c = orasHelmholtz(node, elem, bd, k, parts_c);
    [~, flag, ~, iter] = gmres(A, b, [], gmresTol, maxIter, ap_c);
    its = iter(2) - iter(1);
    if isscalar(nSub), Heff = 1/nSub; else, Heff = 1/max(nSub); end
    fprintf('%-14s %-5s %-5.2f %-6d\n', name, mat2str(nSub), Heff, its);
end

fprintf('\n========== ORAS Parameter Study Complete ==========\n');
