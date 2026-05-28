% VERIFY_LOD_SADDLE_VS_NULLSPACE_SMALL  Diagnostic-only constrained basis check.

fprintf('========== LOD Saddle vs Nullspace Diagnostic ==========\n\n');

fprintf('Test 1: saddle solution equals explicit nullspace solve on tiny system ... ');
rng(1);
n = 8;
nc = 3;
B = randn(n) + 1i * randn(n);
A = B' * B + speye(n);
C = randn(n, nc);
R = randn(n, 2) + 1i * randn(n, 2);

[qSaddle, ~, info] = lodSolveConstrainedSaddle(A, C, R, struct());
Z = null(full(C'));
qBasis = Z * ((Z' * A * Z) \ (Z' * R));

rel = norm(qSaddle - qBasis, 'fro') / max(1, norm(qBasis, 'fro'));
assert(rel < 1e-10, 'Saddle and nullspace diagnostic solves differ: %.3e.', rel);
assert(info.constraintResidual < 1e-10, 'Saddle constraint residual too large.');
fprintf('passed\n');

fprintf('\nNote: explicit nullspace bases are diagnostic-only and are not used by LOD builders.\n');
fprintf('========== LOD saddle/nullspace diagnostic PASSED ==========\n');
