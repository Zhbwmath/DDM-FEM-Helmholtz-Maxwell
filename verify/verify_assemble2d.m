% VERIFY_ASSEMBLE2D  Numerical verification of 2D assembly routines.
%
%   Tests:
%     1. Sparsity pattern & symmetry
%     2. Constant patch test (linear solution reproduced exactly)
%     3. Convergence rate of Poisson solver under uniform refinement

fprintf('========== 2D Assembly Verification ==========\n\n');

%% ---- Test 1: Sparsity & Symmetry ------------------------------------------
fprintf('Test 1: Sparsity pattern and symmetry ... ');

[node, elem] = squaremesh([0, 1, 0, 1], 0.25);
A = assembleStiffness2D(node, elem);
M = assembleMass2D(node, elem);

% Stiffness: symmetric, positive semi-definite (null space = constants)
assert(issymmetric(A), 'Stiffness matrix must be symmetric');
assert(issymmetric(M), 'Mass matrix must be symmetric');

% Check row sums: A * ones = 0 (constant vector in null space)
rowSumA = sum(A, 2);
assert(all(abs(rowSumA) < 1e-12), ...
    'Stiffness matrix must have constant in null space (sum to zero)');

% Mass matrix: all entries positive, diagonal dominance
dA = full(diag(A));
assert(all(dA > 0), 'Stiffness diagonal must be positive');

fprintf('PASSED  (N=%d, NT=%d, nnz(A)=%d)\n', size(node,1), size(elem,1), nnz(A));


%% ---- Test 2: Constant Patch Test ------------------------------------------
fprintf('Test 2: Patch test (constant & linear exact reproduction) ... ');

[nd, el, bd] = squaremesh([0, 1, 0, 1], 0.2);
A = assembleStiffness2D(nd, el);
M = assembleMass2D(nd, el);

% Pure Dirichlet problem: -Delta u = f,  u = u_exact on boundary
% For linear u, f = 0 and stiffness matrix must reproduce u exactly.

u_linear = @(x, y) 1 + 2*x + 3*y;        % Arbitrary linear function
f_linear  = @(x, y) 0;                    % -Delta(linear) = 0

% Assemble RHS:  b = M * f = 0 (or small due to floating point)
b = M * f_linear(nd(:,1), nd(:,2));

% Apply Dirichlet BC for u_linear on all boundary nodes
bdNode = getBoundaryNodes2D(el, bd);
freeNode = setdiff(1:size(nd,1), bdNode)';

u_exact_vals = u_linear(nd(:,1), nd(:,2));
u_bd = u_exact_vals(bdNode);

% Reduced system: A_ff * u_f = b_f - A_fb * u_bd
A_ff = A(freeNode, freeNode);
A_fb = A(freeNode, bdNode);
b_f  = b(freeNode) - A_fb * u_bd;
u_f  = A_ff \ b_f;

% Full solution
uh = zeros(size(nd,1), 1);
uh(bdNode) = u_bd;
uh(freeNode) = u_f;

err = max(abs(uh - u_exact_vals));
assert(err < 1e-10, 'Linear function must be reproduced exactly, error=%.2e', err);

fprintf('PASSED  (max error = %.2e)\n', err);


%% ---- Test 3: Convergence Rate ---------------------------------------------
fprintf('Test 3: Convergence rate (manufactured solution) ...\n');

% Manufactured solution:  u = sin(pi*x)*sin(pi*y)
%   -Delta u = 2*pi^2 * sin(pi*x)*sin(pi*y)
%   u = 0 on boundary of [0,1]^2
u_exact = @(x, y) sin(pi*x) .* sin(pi*y);
f_rhs   = @(x, y) 2*pi^2 * sin(pi*x) .* sin(pi*y);

nRefine = 4;
h_vals    = zeros(nRefine, 1);
errL2     = zeros(nRefine, 1);
errH1     = zeros(nRefine, 1);
dof_vals  = zeros(nRefine, 1);

for k = 1:nRefine
    hk = 2^(-k-1);
    [nd, el, bd] = squaremesh([0, 1, 0, 1], hk);
    Nk = size(nd, 1);

    Ak = assembleStiffness2D(nd, el);
    Mk = assembleMass2D(nd, el);

    % RHS vector
    bx = f_rhs(nd(:,1), nd(:,2));
    bk = Mk * bx;

    % Dirichlet BC (u=0 on boundary)
    bdNodes = getBoundaryNodes2D(el, bd);
    freeNodes = setdiff(1:Nk, bdNodes)';

    b_f = bk(freeNodes);
    A_ff = Ak(freeNodes, freeNodes);
    u_f = A_ff \ b_f;

    uh = zeros(Nk, 1);
    uh(freeNodes) = u_f;

    % Error computation
    u_ex = u_exact(nd(:,1), nd(:,2));
    e_vec = uh - u_ex;

    errL2(k) = sqrt(e_vec' * Mk * e_vec);  % M-weighted L2 norm

    % H1 semi-norm:  |e|_H1^2 = e' * A * e
    errH1(k) = sqrt(e_vec' * Ak * e_vec);

    h_vals(k) = hk;
    dof_vals(k) = Nk;

    if k > 1
        rateL2 = log(errL2(k)/errL2(k-1)) / log(h_vals(k)/h_vals(k-1));
        rateH1 = log(errH1(k)/errH1(k-1)) / log(h_vals(k)/h_vals(k-1));
        fprintf('  h=%.4f  DOF=%5d  |e|_L2=%.4e  rate=%.2f  |e|_H1=%.4e  rate=%.2f\n', ...
            hk, Nk, errL2(k), rateL2, errH1(k), rateH1);
    else
        fprintf('  h=%.4f  DOF=%5d  |e|_L2=%.4e           |e|_H1=%.4e\n', ...
            hk, Nk, errL2(k), errH1(k));
    end
end

% Expected rates: L2 ~ O(h^2), H1 ~ O(h)
assert(rateL2 > 1.80, 'L2 convergence rate %.2f below expected 2.0', rateL2);
assert(rateH1 > 0.80, 'H1 convergence rate %.2f below expected 1.0', rateH1);

fprintf('Test 3: PASSED  (final L2 rate=%.2f, H1 rate=%.2f)\n', rateL2, rateH1);


%% ---- Test 4: Boundary Mass Matrix -----------------------------------------
fprintf('Test 4: Boundary mass matrix ... ');

[nd, el, bd] = squaremesh([0, 1, 0, 1], 0.25);
Mb = assembleBoundaryMass2D(nd, el, bd);

% Should be symmetric, non-negative diagonal
assert(issymmetric(Mb), 'Boundary mass matrix must be symmetric');
dMb = full(diag(Mb));
assert(all(dMb >= 0), 'Boundary mass diagonal must be non-negative');

% Total boundary mass should approximate perimeter length = 4
ones_vec = ones(size(nd,1), 1);
total_bd_mass = ones_vec' * Mb * ones_vec;
assert(abs(total_bd_mass - 4.0) < 0.5, ...
    'Boundary mass integral should approximate perimeter=4, got %.4f', total_bd_mass);

fprintf('PASSED  (total boundary mass = %.4f ~ 4.0)\n', total_bd_mass);

fprintf('\n========== All 2D tests PASSED ==========\n');


% ===========================================================================
function bdNode = getBoundaryNodes2D(elem, bdFlag)
% GETBOUNDARYNODES2D  Return unique list of boundary node indices.
%   A node is a boundary node if it belongs to at least one boundary edge.

bdElem = any(bdFlag == 1, 2);            % elements with at least one bd edge
bdElemIdx = find(bdElem);

bdNode = [];
for k = 1:3                               % check each edge
    % Edges of element t
    switch k
        case 1,  edgeVerts = elem(bdElemIdx, [2, 3]);
        case 2,  edgeVerts = elem(bdElemIdx, [3, 1]);
        case 3,  edgeVerts = elem(bdElemIdx, [1, 2]);
    end
    % Which of these edges are actually boundary edges
    isBd = bdFlag(bdElemIdx, k) == 1;
    bdNode = [bdNode; edgeVerts(isBd, :)]; %#ok<AGROW>
end
bdNode = unique(bdNode(:));
end
