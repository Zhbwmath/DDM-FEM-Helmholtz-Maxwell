% VERIFY_NEDELEC_HELMHOLTZ_DECOMP  Check NE_1 discrete Helmholtz decomposition.

fprintf('========== Nedelec Helmholtz Decomposition ==========\n\n');

tol = 5e-11;
semiTol = 1e-7;

% 2D: exact gradients must be recovered and have zero discrete curl.
[node2, elem2] = squaremesh([0, 1, 0, 1], 1/4);
[G2, ~] = nedelecGradientMatrix(node2, elem2);
M2 = assembleNedMass2D(node2, elem2);
C2 = assembleCurlCurl2D(node2, elem2);
opts2 = struct('massMatrix', M2, 'curlCurlMatrix', C2);
phi2 = sin(node2(:,1) + 2*node2(:,2));
u2 = G2 * phi2;

[uGrad2, uComp2, ~, info2] = nedelecHelmholtzDecomp(node2, elem2, u2, ...
    struct('massMatrix', M2));
[curlGrad2, ~] = nedelecCurlSeminorm(node2, elem2, uGrad2, opts2);
[divComp2, divInfo2] = nedelecDivSeminorm(node2, elem2, uComp2, opts2);
assertRelative('2D gradient recovery', uGrad2, u2, tol);
assert(norm(uComp2) <= tol * max(norm(u2), 1), '2D gradient complement is not zero.');
assert(curlGrad2 <= semiTol * max(norm(uGrad2), 1), '2D gradient component is not curl-free.');
assert(divComp2 <= semiTol * max(norm(u2), 1), '2D complement is not weakly divergence-free.');
assert(info2.reconstructionError <= eps, '2D reconstruction is not exact.');

% 2D: arbitrary fields split into a gradient and an M-orthogonal complement.
uRand2 = sin((1:size(G2,1))' * 0.37) + cos((1:size(G2,1))' * 0.11);
[uGrad2b, uComp2b, ~, info2b] = nedelecHelmholtzDecomp(node2, elem2, uRand2, ...
    struct('massMatrix', M2));
[curlGrad2b, ~] = nedelecCurlSeminorm(node2, elem2, uGrad2b, opts2);
[divComp2b, divInfo2b] = nedelecDivSeminorm(node2, elem2, uComp2b, opts2);
assertRelative('2D reconstruction', uGrad2b + uComp2b, uRand2, eps);
assert(info2b.orthogonalityNorm <= tol * max(norm(G2' * (M2 * uRand2)), 1), ...
    '2D complement is not M-orthogonal to gradients.');
assert(curlGrad2b <= semiTol * max(norm(uGrad2b), 1), '2D projected gradient is not curl-free.');
assert(divComp2b <= semiTol * max(norm(uRand2), 1), ...
    '2D projected complement is not weakly divergence-free.');

% 3D: same exact-gradient and orthogonality checks on a tetrahedral mesh.
[node3, elem3] = cubemesh([0, 1, 0, 1, 0, 1], 1/2);
[G3, ~] = nedelecGradientMatrix(node3, elem3);
M3 = assembleNedMass3D(node3, elem3);
C3 = assembleCurlCurl3D(node3, elem3);
opts3 = struct('massMatrix', M3, 'curlCurlMatrix', C3);
phi3 = node3(:,1) - 2*node3(:,2) + 0.5*sin(node3(:,3));
u3 = G3 * phi3;

[uGrad3, uComp3, ~, info3] = nedelecHelmholtzDecomp(node3, elem3, u3, ...
    struct('massMatrix', M3));
[curlGrad3, ~] = nedelecCurlSeminorm(node3, elem3, uGrad3, opts3);
[divComp3, divInfo3] = nedelecDivSeminorm(node3, elem3, uComp3, opts3);
assertRelative('3D gradient recovery', uGrad3, u3, tol);
assert(norm(uComp3) <= tol * max(norm(u3), 1), '3D gradient complement is not zero.');
assert(curlGrad3 <= semiTol * max(norm(uGrad3), 1), '3D gradient component is not curl-free.');
assert(divComp3 <= semiTol * max(norm(u3), 1), '3D complement is not weakly divergence-free.');
assert(info3.reconstructionError <= eps, '3D reconstruction is not exact.');

uRand3 = sin((1:size(G3,1))' * 0.23) - cos((1:size(G3,1))' * 0.41);
[uGrad3b, uComp3b, ~, info3b] = nedelecHelmholtzDecomp(node3, elem3, uRand3, ...
    struct('massMatrix', M3));
[curlGrad3b, ~] = nedelecCurlSeminorm(node3, elem3, uGrad3b, opts3);
[divComp3b, divInfo3b] = nedelecDivSeminorm(node3, elem3, uComp3b, opts3);
assertRelative('3D reconstruction', uGrad3b + uComp3b, uRand3, eps);
assert(info3b.orthogonalityNorm <= tol * max(norm(G3' * (M3 * uRand3)), 1), ...
    '3D complement is not M-orthogonal to gradients.');
assert(curlGrad3b <= semiTol * max(norm(uGrad3b), 1), '3D projected gradient is not curl-free.');
assert(divComp3b <= semiTol * max(norm(uRand3), 1), ...
    '3D projected complement is not weakly divergence-free.');

fprintf('2D curl(uGrad) exact-gradient seminorm: %.3e\n', curlGrad2);
fprintf('2D div_h(uComp) exact-gradient seminorm: %.3e (raw %.3e)\n', ...
    divComp2, divInfo2.weakDivergenceResidualNorm);
fprintf('2D curl(uGrad) projected seminorm: %.3e\n', curlGrad2b);
fprintf('2D div_h(uComp) projected seminorm: %.3e (raw %.3e)\n', ...
    divComp2b, divInfo2b.weakDivergenceResidualNorm);
fprintf('2D orthogonality: %.3e\n', info2b.orthogonalityNorm);
fprintf('3D curl(uGrad) exact-gradient seminorm: %.3e\n', curlGrad3);
fprintf('3D div_h(uComp) exact-gradient seminorm: %.3e (raw %.3e)\n', ...
    divComp3, divInfo3.weakDivergenceResidualNorm);
fprintf('3D curl(uGrad) projected seminorm: %.3e\n', curlGrad3b);
fprintf('3D div_h(uComp) projected seminorm: %.3e (raw %.3e)\n', ...
    divComp3b, divInfo3b.weakDivergenceResidualNorm);
fprintf('3D orthogonality: %.3e\n', info3b.orthogonalityNorm);
fprintf('========== Done ==========\n');

function assertRelative(name, got, expected, tol)
err = norm(got - expected, 'fro') / max(norm(expected, 'fro'), 1);
assert(err <= tol, '%s failed with relative error %.3e.', name, err);
end
