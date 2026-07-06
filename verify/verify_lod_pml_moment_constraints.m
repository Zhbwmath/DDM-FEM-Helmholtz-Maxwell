% VERIFY_LOD_PML_MOMENT_CONSTRAINTS  Checks exact L2 moment constraints.

fprintf('========== LOD PML Moment Constraint Verification ==========\n\n');

[nodeH, elemH, ~] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);
P = prolongateNestedP1(nodeH, elemH, nodeh);
M = assembleMass2D(nodeh, elemh, 1);
Crows = P' * M;
MH = P' * M * P;

fprintf('Test 1: coarse mass consistency ... ');
relMH = norm(MH - P' * M * P, 'fro') / max(1, norm(MH, 'fro'));
assert(relMH < 1e-14, 'Coarse mass identity failed.');
fprintf('PASSED\n');

fprintf('Test 2: projected random vector is in ker(P''M) ... ');
rng(12);
w = randn(size(nodeh, 1), 1);
PiW = P * (MH \ (P' * M * w));
wk = w - PiW;
kerResidual = norm(Crows * wk) / max(1, norm(wk));
projResidual = sqrt(real(PiW' * M * PiW)) / max(1, sqrt(real(w' * M * w)));
assert(kerResidual < 1e-12, 'Moment-kernel residual too large %.3e.', kerResidual);
assert(projResidual > 1e-8, 'Projection check did not exercise a nonzero coarse part.');
fprintf('PASSED (ker residual %.2e)\n', kerResidual);

fprintf('Test 3: local constraint assembly matches restricted global rows ... ');
patch = lodBuildPatches(nodeH, elemH, nodeh, elemh, bdh, 1);
opts = struct('constraintTolerance', 1e-12);
T = 1;
[C, info] = lodMomentConstraints(Crows, patch, T, opts, [], []);
sub = lodGetPatchSubmesh(patch, T);
freeGlobal = sub.local2global(sub.freeLocalDof);
Cref = Crows(info.coarseDof, freeGlobal)';
relC = norm(C - Cref, 'fro') / max(1, norm(Cref, 'fro'));
assert(relC < 1e-14, 'Patch moment constraints mismatch %.3e.', relC);
fprintf('PASSED (rel %.2e)\n', relC);

fprintf('\n========== LOD PML moment constraint tests PASSED ==========\n');
