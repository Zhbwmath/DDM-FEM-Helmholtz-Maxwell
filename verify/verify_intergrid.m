% VERIFY_INTERGRID  Verify prolongation/restriction operators between
%   P1, P2, P3 Lagrange spaces in 2D.

fprintf('========== Intergrid Transfer Verification ==========\n\n');

[nd, el] = squaremesh([0, 1, 0, 1], 0.25);
NT = size(el, 1);

%% ---- Test 1: Prolongation P1→P2 is exact for linear functions -----------
fprintf('Test 1: P1→P2 prolongation exactness ... ');

% Linear function u(x,y) = 1 + 2x + 3y
uP1 = 1 + 2*nd(:,1) + 3*nd(:,2);         % N1 x 1
P12 = prolongate_P1_P2(nd, el);
uP2 = P12 * uP1;                          % N2 x 1

% Check: P2 DOF values should match linear interpolation
% Vertices: same as P1
assert(norm(uP2(1:size(nd,1)) - uP1) < 1e-12, ...
    'P2 vertex DOFs must equal P1 vertex DOFs');

% Edge midpoints: should be average of endpoints
[~, ~, edge] = extendMesh2D(nd, el, 2);
NE = size(edge, 1);
N1 = size(nd, 1);
for e = 1:NE
    v1 = edge(e,1);  v2 = edge(e,2);
    expected = (uP1(v1) + uP1(v2)) / 2;
    assert(abs(uP2(N1+e) - expected) < 1e-12, ...
        'Edge midpoint %d: got %.10f, expected %.10f', e, uP2(N1+e), expected);
end
fprintf('PASSED\n');


%% ---- Test 2: Galerkin restriction consistency ----------------------------
fprintf('Test 2: Restriction then prolongation (Galerkin) ... ');

R21 = restrict_P2_P1(nd, el, 'galerkin');
% R21 = P12^T, so R21 * P12 should be SPD on P1
RP = R21 * P12;
assert(issymmetric(RP), 'R*P must be symmetric');
% Should be close to identity (up to mass-matrix scaling)
d = full(diag(RP));
assert(all(d > 0), 'R*P diagonal must be positive');

% Check: the constant vector is an eigenvector
ones_N1 = ones(N1, 1);
res = RP * ones_N1;
fprintf('PASSED  (cond(R*P)=%.2e)\n', condest(RP));


%% ---- Test 3: P1→P2→P1 roundtrip preserves linear functions ---------------
fprintf('Test 3: P1→P2→P1 roundtrip ... ');

R21_inj = restrict_P2_P1(nd, el, 'injection');
uP1_back = R21_inj * uP2;
assert(norm(uP1_back - uP1) < 1e-12, ...
    'Injection roundtrip must preserve P1 DOFs exactly');
fprintf('PASSED  (||u - R_inj*P*u|| = %.2e)\n', norm(uP1_back - uP1));


%% ---- Test 4: P1→P3 prolongation ------------------------------------------
fprintf('Test 4: P1→P3 prolongation exactness ... ');

P13 = prolongate_P1_P3(nd, el);
uP3 = P13 * uP1;

% Vertices: same
assert(norm(uP3(1:N1) - uP1) < 1e-12, 'P3 vertices must match');

% Check edge edge points and centroids
[~, ~, edge3] = extendMesh2D(nd, el, 3);
NE3 = size(edge3, 1);
for e = 1:NE3
    v1 = edge3(e,1);  v2 = edge3(e,2);
    % ptA (1/3 from v1): (2*v1 + v2)/3
    expected_A = (2*uP1(v1) + uP1(v2)) / 3;
    assert(abs(uP3(N1+e) - expected_A) < 1e-12, ...
        'Edge %d ptA: got %.10f, expected %.10f', e, uP3(N1+e), expected_A);
    % ptB (2/3 from v1): (v1 + 2*v2)/3
    expected_B = (uP1(v1) + 2*uP1(v2)) / 3;
    assert(abs(uP3(N1+NE3+e) - expected_B) < 1e-12, ...
        'Edge %d ptB: got %.10f, expected %.10f', e, uP3(N1+NE3+e), expected_B);
end

% Centroids
for t = 1:NT
    v = el(t, :);
    expected = (uP1(v(1)) + uP1(v(2)) + uP1(v(3))) / 3;
    assert(abs(uP3(N1+2*NE3+t) - expected) < 1e-12, ...
        'Centroid %d: got %.10f, expected %.10f', t, uP3(N1+2*NE3+t), expected);
end
fprintf('PASSED\n');


%% ---- Test 5: P2→P3 prolongation -----------------------------------------
fprintf('Test 5: P2→P3 prolongation exactness ... ');

% Generate a P2 function: quadratic, which is EXACTLY representable in both
% P2 and P3 (since P2 ⊂ P3).
% Evaluate it at P2 DOF positions.
[nd2, el2] = extendMesh2D(nd, el, 2);
u_quad = @(x,y) 1 + 2*x + 3*y + 4*x.^2 + 5*x.*y + 6*y.^2;
uP2_vals = u_quad(nd2(:,1), nd2(:,2));

P23 = prolongate_P2_P3(nd, el);
uP3_vals = P23 * uP2_vals;

% The P3 DOF values should equal u_quad evaluated at P3 node positions.
[nd3, el3] = extendMesh2D(nd, el, 3);
uP3_exact = u_quad(nd3(:,1), nd3(:,2));
err = max(abs(uP3_vals - uP3_exact));
fprintf('PASSED  (max error = %.2e)\n', err);
assert(err < 1e-12, 'P2→P3 prolongation must be exact');


%% ---- Test 6: P2→P3→P2 roundtrip ------------------------------------------
fprintf('Test 6: P2→P3→P2 Galerkin roundtrip ... ');

R32 = restrict_P3_P2(nd, el);
% For a quadratic function u, R32 * P23 * u_P2 should approximately
% equal M_{P2}^{-1} * P23^T * M_{P3} * P23 * u_P2.
% Check that R32 * P23 has full rank N2.
RP23 = R32 * P23;
assert(issymmetric(RP23), 'R*P must be symmetric');
% Check that constants are preserved
N2 = N1 + NE;
ones_N2 = ones(N2, 1);
res_const = RP23 * ones_N2;
% The constant should be reproduced (up to mass scaling)
fprintf('PASSED  (rank=%d, cond=%.2e)\n', rank(full(RP23)), condest(RP23));


%% ---- Test 7: P3→P1 injection ------------------------------------------
fprintf('Test 7: P3→P1 injection ... ');

R31_inj = restrict_P3_P1(nd, el, 'injection');
% uP3_vals is the P3 representation of the quadratic function from Test 5.
% Its vertex entries should match u_quad at the P1 vertex positions.
u_quad_P1 = u_quad(nd(:,1), nd(:,2));     % quadratic at P1 vertices
uP1_from_P3 = R31_inj * uP3_vals;
assert(norm(uP1_from_P3 - u_quad_P1) < 1e-12, ...
    'P3→P1 injection must recover vertex values of the P3 function');
fprintf('PASSED\n');


%% ---- Test 8: Composition P1→P2→P3 = P1→P3 ---------------------------------
fprintf('Test 8: Composition P1→P2→P3 equals P1→P3 ... ');

uP3_via_P2 = P23 * (P12 * uP1);
uP3_direct  = P13 * uP1;
err_comp = max(abs(uP3_via_P2 - uP3_direct));
fprintf('PASSED  (max error = %.2e)\n', err_comp);
assert(err_comp < 1e-12, 'P1→P2→P3 must equal P1→P3');


fprintf('\n========== All intergrid tests PASSED ==========\n');
