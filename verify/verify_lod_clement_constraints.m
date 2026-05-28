% VERIFY_LOD_CLEMENT_CONSTRAINTS  Check Clement patch constraints for LOD.

fprintf('========== LOD Clement Constraint Verification ==========\n\n');

[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);
patch = lodBuildPatches2D(nodeH, elemH, nodeh, elemh, bdh, 1);
Q = weightedClementP1(nodeh, elemh, nodeH, elemH);

fprintf('Test 1: local constraint dimensions and active rows ... ');
for T = 1:size(elemH, 1)
    [C, info] = lodClementConstraints(Q, patch, T, struct());
    nFree = nnz(patch.freeLocalDof{T});
    assert(size(C, 1) == nFree, 'Constraint row count must match free local DOFs.');
    assert(size(C, 2) == numel(info.coarseDof), 'Constraint info mismatch.');
    assert(~isempty(info.coarseDof), 'Each patch should have active Clement rows.');
end
fprintf('passed\n');

fprintf('Test 2: saddle solution satisfies Clement constraints ... ');
T = 1;
[C, ~] = lodClementConstraints(Q, patch, T, struct());
n = size(C, 1);
A = speye(n) + 0.1 * sprandsym(n, 0.4, 0.1, 1);
R = ones(n, 2);
[q, ~, info] = lodSolveConstrainedSaddle(A, C, R, struct());
assert(info.constraintResidual < 1e-10, 'Saddle solve did not satisfy C''q=0.');
assert(norm(C' * q, 'fro') < 1e-10 * max(1, norm(q, 'fro')), 'Constraint residual too large.');
fprintf('passed\n');

fprintf('\n========== LOD Clement constraint tests PASSED ==========\n');
