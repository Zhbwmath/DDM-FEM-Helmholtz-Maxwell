% VERIFY_NEDELEC_HELMHOLTZ_DECOMP  Check NE_1 discrete Helmholtz decomposition.

fprintf('========== Nedelec Helmholtz Decomposition ==========\n\n');

tol = 5e-11;

% 2D: exact gradients must be recovered and have zero discrete curl.
[node2, elem2] = squaremesh([0, 1, 0, 1], 1/4);
[G2, ~] = nedelecGradientMatrix(node2, elem2);
M2 = assembleNedMass2D(node2, elem2);
C2 = assembleCurlCurl2D(node2, elem2);
phi2 = sin(node2(:,1) + 2*node2(:,2));
u2 = G2 * phi2;

[uGrad2, uComp2, ~, info2] = nedelecHelmholtzDecomp(node2, elem2, u2, ...
    struct('massMatrix', M2));
assertRelative('2D gradient recovery', uGrad2, u2, tol);
assert(norm(uComp2) <= tol * max(norm(u2), 1), '2D gradient complement is not zero.');
assert(norm(C2 * uGrad2, inf) <= tol * max(norm(uGrad2), 1), ...
    '2D gradient component is not curl-free.');
assert(info2.reconstructionError <= eps, '2D reconstruction is not exact.');

% 2D: arbitrary fields split into a gradient and an M-orthogonal complement.
uRand2 = sin((1:size(G2,1))' * 0.37) + cos((1:size(G2,1))' * 0.11);
[uGrad2b, uComp2b, ~, info2b] = nedelecHelmholtzDecomp(node2, elem2, uRand2, ...
    struct('massMatrix', M2));
assertRelative('2D reconstruction', uGrad2b + uComp2b, uRand2, eps);
assert(info2b.orthogonalityNorm <= tol * max(norm(G2' * (M2 * uRand2)), 1), ...
    '2D complement is not M-orthogonal to gradients.');

% 3D: same exact-gradient and orthogonality checks on a tetrahedral mesh.
[node3, elem3] = cubemesh([0, 1, 0, 1, 0, 1], 1/2);
[G3, ~] = nedelecGradientMatrix(node3, elem3);
M3 = assembleNedMass3D(node3, elem3);
C3 = assembleCurlCurl3D(node3, elem3);
phi3 = node3(:,1) - 2*node3(:,2) + 0.5*sin(node3(:,3));
u3 = G3 * phi3;

[uGrad3, uComp3, ~, info3] = nedelecHelmholtzDecomp(node3, elem3, u3, ...
    struct('massMatrix', M3));
assertRelative('3D gradient recovery', uGrad3, u3, tol);
assert(norm(uComp3) <= tol * max(norm(u3), 1), '3D gradient complement is not zero.');
assert(norm(C3 * uGrad3, inf) <= tol * max(norm(uGrad3), 1), ...
    '3D gradient component is not curl-free.');
assert(info3.reconstructionError <= eps, '3D reconstruction is not exact.');

uRand3 = sin((1:size(G3,1))' * 0.23) - cos((1:size(G3,1))' * 0.41);
[uGrad3b, uComp3b, ~, info3b] = nedelecHelmholtzDecomp(node3, elem3, uRand3, ...
    struct('massMatrix', M3));
assertRelative('3D reconstruction', uGrad3b + uComp3b, uRand3, eps);
assert(info3b.orthogonalityNorm <= tol * max(norm(G3' * (M3 * uRand3)), 1), ...
    '3D complement is not M-orthogonal to gradients.');

fprintf('2D orthogonality: %.3e\n', info2b.orthogonalityNorm);
fprintf('3D orthogonality: %.3e\n', info3b.orthogonalityNorm);
fprintf('========== Done ==========\n');

function assertRelative(name, got, expected, tol)
err = norm(got - expected, 'fro') / max(norm(expected, 'fro'), 1);
assert(err <= tol, '%s failed with relative error %.3e.', name, err);
end
