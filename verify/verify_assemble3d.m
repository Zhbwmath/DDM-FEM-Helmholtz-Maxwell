% VERIFY_ASSEMBLE3D  Numerical verification of 3D assembly routines.
%
%   Tests:
%     1. Sparsity pattern & symmetry
%     2. Constant patch test (linear solution reproduced exactly)
%     3. Convergence rate of Poisson solver under uniform refinement

fprintf('========== 3D Assembly Verification ==========\n\n');

%% ---- Test 1: Sparsity & Symmetry ------------------------------------------
fprintf('Test 1: Sparsity pattern and symmetry ... ');

[node, elem] = cubemesh([0, 1, 0, 1, 0, 1], 0.4);
A = assembleStiffness3D(node, elem);
M = assembleMass3D(node, elem);

assert(issymmetric(A), 'Stiffness matrix must be symmetric');
assert(issymmetric(M), 'Mass matrix must be symmetric');

% Constant vector in null space of stiffness
rowSumA = sum(A, 2);
assert(all(abs(rowSumA) < 1e-12), ...
    'Stiffness matrix must have constant in null space');

dA = full(diag(A));
assert(all(dA > 0), 'Stiffness diagonal must be positive');

fprintf('PASSED  (N=%d, NT=%d, nnz(A)=%d)\n', size(node,1), size(elem,1), nnz(A));


%% ---- Test 2: Constant Patch Test ------------------------------------------
fprintf('Test 2: Patch test (constant & linear exact reproduction) ... ');

[nd, el, bd] = cubemesh([0, 1, 0, 1, 0, 1], 0.4);
A = assembleStiffness3D(nd, el);
M = assembleMass3D(nd, el);

% Linear function:  -Delta u = 0, should be reproduced exactly (up to BC)
u_linear = @(x, y, z) 1 + 2*x + 3*y + 4*z;

b = M * zeros(size(nd,1), 1);     % f = 0

% Dirichlet BC on all boundary nodes
bdNode = getBoundaryNodes3D(el, bd);
freeNode = setdiff(1:size(nd,1), bdNode)';

u_exact_vals = u_linear(nd(:,1), nd(:,2), nd(:,3));
u_bd = u_exact_vals(bdNode);

% Reduced solve
u_f = A(freeNode, freeNode) \ (b(freeNode) - A(freeNode, bdNode) * u_bd);
uh = zeros(size(nd,1), 1);
uh(bdNode) = u_bd;
uh(freeNode) = u_f;

err = max(abs(uh - u_exact_vals));
assert(err < 1e-10, 'Linear function must be reproduced exactly, error=%.2e', err);

fprintf('PASSED  (max error = %.2e)\n', err);


%% ---- Test 3: Convergence Rate ---------------------------------------------
fprintf('Test 3: Convergence rate (manufactured solution) ...\n');

% Manufactured solution:  u = sin(pi*x)*sin(pi*y)*sin(pi*z)
%   -Delta u = 3*pi^2 * sin(pi*x)*sin(pi*y)*sin(pi*z)
%   u = 0 on boundary of [0,1]^3
u_exact = @(x, y, z) sin(pi*x) .* sin(pi*y) .* sin(pi*z);
f_rhs   = @(x, y, z) 3*pi^2 * sin(pi*x) .* sin(pi*y) .* sin(pi*z);

nRefine = 3;                              % 3D is heavy; 3 refinements enough
h_vals    = zeros(nRefine, 1);
errL2     = zeros(nRefine, 1);
errH1     = zeros(nRefine, 1);
dof_vals  = zeros(nRefine, 1);

for k = 1:nRefine
    hk = 2^(-k-1);                        % h = 1/4, 1/8, 1/16
    [nd, el, bd] = cubemesh([0, 1, 0, 1, 0, 1], hk);
    Nk = size(nd, 1);

    fprintf('  Assembling h=%.4f (N=%d, NT=%d) ... ', hk, Nk, size(el,1));

    Ak = assembleStiffness3D(nd, el);
    Mk = assembleMass3D(nd, el);

    % RHS
    bx = f_rhs(nd(:,1), nd(:,2), nd(:,3));
    bk = Mk * bx;

    % Dirichlet BC
    bdNodes = getBoundaryNodes3D(el, bd);
    freeNodes = setdiff(1:Nk, bdNodes)';

    u_f = Ak(freeNodes, freeNodes) \ bk(freeNodes);
    uh = zeros(Nk, 1);
    uh(freeNodes) = u_f;

    u_ex = u_exact(nd(:,1), nd(:,2), nd(:,3));
    e_vec = uh - u_ex;

    errL2(k) = sqrt(e_vec' * Mk * e_vec);
    errH1(k) = sqrt(e_vec' * Ak * e_vec);
    h_vals(k) = hk;
    dof_vals(k) = Nk;

    if k > 1
        rateL2 = log(errL2(k)/errL2(k-1)) / log(h_vals(k)/h_vals(k-1));
        rateH1 = log(errH1(k)/errH1(k-1)) / log(h_vals(k)/h_vals(k-1));
        fprintf('|e|_L2=%.4e  rate=%.2f  |e|_H1=%.4e  rate=%.2f\n', ...
            errL2(k), rateL2, errH1(k), rateH1);
    else
        fprintf('|e|_L2=%.4e           |e|_H1=%.4e\n', errL2(k), errH1(k));
    end
end

assert(rateL2 > 1.80, 'L2 convergence rate %.2f below expected 2.0', rateL2);
assert(rateH1 > 0.80, 'H1 convergence rate %.2f below expected 1.0', rateH1);

fprintf('Test 3: PASSED  (final L2 rate=%.2f, H1 rate=%.2f)\n', rateL2, rateH1);


%% ---- Test 4: Boundary Mass Matrix -----------------------------------------
fprintf('Test 4: Boundary mass matrix ... ');

[nd, el, bd] = cubemesh([0, 1, 0, 1, 0, 1], 0.4);
Mb = assembleBoundaryMass3D(nd, el, bd);

assert(issymmetric(Mb), 'Boundary mass matrix must be symmetric');
dMb = full(diag(Mb));
assert(all(dMb >= 0), 'Boundary mass diagonal must be non-negative');

% Total boundary mass should approximate surface area = 6
ones_vec = ones(size(nd,1), 1);
total_bd_mass = ones_vec' * Mb * ones_vec;
assert(abs(total_bd_mass - 6.0) < 1.0, ...
    'Boundary mass integral should approximate surface area=6, got %.4f', total_bd_mass);

fprintf('PASSED  (total boundary mass = %.4f ~ 6.0)\n', total_bd_mass);

fprintf('\n========== All 3D tests PASSED ==========\n');


% ===========================================================================
function bdNode = getBoundaryNodes3D(elem, bdFlag)
% GETBOUNDARYNODES3D  Return unique list of boundary node indices.
%   A node is on the Dirichlet boundary if it belongs to at least one
%   boundary face.

faceVerts = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
bdNode = [];

for f = 1:4
    isFace = bdFlag(:,f) == 1;
    if ~any(isFace), continue; end
    fv = faceVerts{f};
    bdNode = [bdNode; elem(isFace, fv(1)); elem(isFace, fv(2)); elem(isFace, fv(3))]; %#ok<AGROW>
end
bdNode = unique(bdNode);
end
