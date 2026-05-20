% VERIFY_ASSEMBLE1D  Numerical verification of 1D assembly routines.
%
%   Tests:
%     1. Sparsity pattern & symmetry
%     2. Constant patch test (linear solution reproduced exactly)
%     3. Convergence rate of Poisson solver under uniform refinement
%     4. Boundary mass matrix
%     5. Higher-order convergence (P2, P3)

fprintf('========== 1D Assembly Verification ==========\n\n');

%% ---- Test 1: Sparsity & Symmetry ------------------------------------------
fprintf('Test 1: Sparsity pattern and symmetry ... ');

[node, elem] = linemesh(0, 1, 16);
A = assembleStiffness1D(node, elem);
M = assembleMass1D(node, elem);

assert(issymmetric(A), 'Stiffness matrix must be symmetric');
assert(issymmetric(M), 'Mass matrix must be symmetric');

% Stiffness: constant in null space (rigid body translation)
rowSumA = full(sum(A, 2));
assert(all(abs(rowSumA(2:end-1)) < 1e-12), ...
    'Interior rows of stiffness must sum to zero');

dA = full(diag(A));
assert(all(dA > 0), 'Stiffness diagonal must be positive');

fprintf('PASSED  (N=%d, NT=%d, nnz(A)=%d)\n', size(node,1), size(elem,1), nnz(A));


%% ---- Test 2: Patch Test (linear exact reproduction) -----------------------
fprintf('Test 2: Patch test (linear exact reproduction) ... ');

[node, elem, bd] = linemesh(0, 1, 10);
A = assembleStiffness1D(node, elem);
M = assembleMass1D(node, elem);

% Linear function: u(x) = 2 + 3x,  -u'' = 0
u_linear = @(x) 2 + 3*x;
f_linear  = @(x) 0 * x;

b = M * f_linear(node);

% Dirichlet BC on both endpoints (map local vertex flags to global nodes)
bdNodes = [elem(bd(:,1)==1, 1); elem(bd(:,2)==1, 2)];
bdNodes = unique(bdNodes);
freeNodes = setdiff(1:size(node,1), bdNodes)';

u_ex = u_linear(node);
u_bd = u_ex(bdNodes);

A_ff = A(freeNodes, freeNodes);
A_fb = A(freeNodes, bdNodes);
b_f  = b(freeNodes) - A_fb * u_bd;
u_f  = A_ff \ b_f;

uh = zeros(size(node,1), 1);
uh(bdNodes) = u_bd;
uh(freeNodes) = u_f;

err = max(abs(uh - u_ex));
assert(err < 1e-10, 'Linear function must be reproduced exactly, error=%.2e', err);

fprintf('PASSED  (max error = %.2e)\n', err);


%% ---- Test 3: Convergence Rate ---------------------------------------------
fprintf('Test 3: Convergence rate (manufactured solution) ...\n');

% u = sin(pi*x),  -u'' = pi^2 * sin(pi*x),  u(0)=u(1)=0
u_exact = @(x) sin(pi*x);
f_rhs   = @(x) pi^2 * sin(pi*x);

nRefine = 5;
hVals    = zeros(nRefine, 1);
errL2    = zeros(nRefine, 1);
errH1    = zeros(nRefine, 1);
dofVals  = zeros(nRefine, 1);

for k = 1:nRefine
    nElem = 2^(k+2);
    [nd, el, bd] = linemesh(0, 1, nElem);
    Nk = size(nd, 1);

    Ak = assembleStiffness1D(nd, el);
    Mk = assembleMass1D(nd, el);

    bk = Mk * f_rhs(nd);

    bdNodes = [el(bd(:,1)==1, 1); el(bd(:,2)==1, 2)];
    bdNodes = unique(bdNodes);
    freeNodes = setdiff(1:Nk, bdNodes)';

    b_f = bk(freeNodes);
    A_ff = Ak(freeNodes, freeNodes);
    u_f = A_ff \ b_f;

    uh = zeros(Nk, 1);
    uh(freeNodes) = u_f;

    u_ex = u_exact(nd);
    e_vec = uh - u_ex;

    errL2(k) = sqrt(e_vec' * Mk * e_vec);
    errH1(k) = sqrt(e_vec' * Ak * e_vec);

    hVals(k) = 1 / nElem;
    dofVals(k) = Nk;

    if k > 1
        rateL2 = log(errL2(k)/errL2(k-1)) / log(hVals(k)/hVals(k-1));
        rateH1 = log(errH1(k)/errH1(k-1)) / log(hVals(k)/hVals(k-1));
        fprintf('  h=%.4f  DOF=%5d  |e|_L2=%.4e  rate=%.2f  |e|_H1=%.4e  rate=%.2f\n', ...
            hVals(k), Nk, errL2(k), rateL2, errH1(k), rateH1);
    else
        fprintf('  h=%.4f  DOF=%5d  |e|_L2=%.4e           |e|_H1=%.4e\n', ...
            hVals(k), Nk, errL2(k), errH1(k));
    end
end

assert(rateL2 > 1.80, 'L2 convergence rate %.2f below expected 2.0', rateL2);
assert(rateH1 > 0.80, 'H1 convergence rate %.2f below expected 1.0', rateH1);

fprintf('Test 3: PASSED  (final L2 rate=%.2f, H1 rate=%.2f)\n', rateL2, rateH1);


%% ---- Test 4: Boundary Mass Matrix -----------------------------------------
fprintf('Test 4: Boundary mass matrix ... ');

[nd, el, bd] = linemesh(0, 1, 16);
Mb = assembleBoundaryMass1D(nd, el, bd);

assert(issymmetric(Mb), 'Boundary mass matrix must be symmetric');
dMb = full(diag(Mb));
assert(all(dMb >= 0), 'Boundary mass diagonal must be non-negative');

% For P1: Mb has 1 at endpoints (nodes 1 and N), 0 elsewhere
assert(abs(dMb(1) - 1) < 1e-12 && abs(dMb(end) - 1) < 1e-12, ...
    'Boundary mass diagonal should be 1 at endpoints for P1');

% Sum should be 2 (two endpoints, each weight 1)
total_bd = full(sum(Mb, 'all'));
assert(abs(total_bd - 2) < 1e-12, ...
    'Total boundary mass should be 2, got %.4f', total_bd);

fprintf('PASSED  (total boundary mass = %.4f = 2.0)\n', total_bd);


%% ---- Test 5: Higher-order convergence (P2, P3) ----------------------------
fprintf('Test 5: Higher-order convergence ...\n');

for deg = [2, 3]
    errL2_ho = zeros(3, 1);
    hVals_ho = zeros(3, 1);
    for k = 1:3
        nElem = 2^(k+1);
        [nd, el, bd] = linemesh(0, 1, nElem);

        Ak = assembleStiffness1D(nd, el, deg);
        Mk = assembleMass1D(nd, el, deg);

        Nk = size(Ak, 1);

        % Evaluate f at all nodes (mesh may have been extended)
        bk = Mk * f_rhs(linspace(0, 1, Nk)');

        % Boundary nodes: x=0 and x=1
        node_coords = linspace(0, 1, Nk)';
        bdNodes = find(abs(node_coords - 0) < 1e-12 | abs(node_coords - 1) < 1e-12);
        freeNodes = setdiff(1:Nk, bdNodes)';

        b_f = bk(freeNodes);
        A_ff = Ak(freeNodes, freeNodes);
        u_f = A_ff \ b_f;

        uh = zeros(Nk, 1);
        uh(freeNodes) = u_f;

        u_ex = u_exact(node_coords);
        e_vec = uh - u_ex;

        errL2_ho(k) = sqrt(e_vec' * Mk * e_vec);
        hVals_ho(k) = 1 / nElem;

        if k > 1
            rate = log(errL2_ho(k)/errL2_ho(k-1)) / log(hVals_ho(k)/hVals_ho(k-1));
            fprintf('  P%d: h=%.4f  DOF=%5d  |e|_L2=%.4e  rate=%.2f\n', ...
                deg, hVals_ho(k), Nk, errL2_ho(k), rate);
        else
            fprintf('  P%d: h=%.4f  DOF=%5d  |e|_L2=%.4e\n', ...
                deg, hVals_ho(k), Nk, errL2_ho(k));
        end
    end
    expectedRate = deg + 1;
    assert(rate > expectedRate - 0.30, ...
        'P%d L2 convergence rate %.2f below expected %d', deg, rate, expectedRate);
end

fprintf('Test 5: PASSED\n');

fprintf('\n========== All 1D tests PASSED =========\n');
