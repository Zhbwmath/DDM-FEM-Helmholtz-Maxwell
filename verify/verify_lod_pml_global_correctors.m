% VERIFY_LOD_PML_GLOBAL_CORRECTORS  Checks global L2-moment PML correctors.

fprintf('========== LOD PML Global Corrector Verification ==========\n\n');

[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);
k = 3;
pml = struct('physicalBox', [0, 1, 0, 1], 'pmlBox', [0, 1, 0, 1], ...
    'sigmaMax', 0, 'quadOrder', 4);
opts = struct('constraintTolerance', 1e-12, 'solverMode', 'direct', ...
    'degree', 1);
problem = helmholtzPMLLODProblem2D(nodeH, elemH, bdH, nodeh, elemh, bdh, ...
    k, pml, 1, opts);
basisId = 5;

fprintf('Test 1: global primal corrector variational orthogonality ... ');
opts.basisSide = 'trial';
trial = lodCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisId, Inf, opts);
checkOrthogonality(problem, trial, basisId, 'trial');
fprintf('PASSED\n');

fprintf('Test 2: global adjoint corrector variational orthogonality ... ');
opts.basisSide = 'test';
test = lodCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisId, Inf, opts);
checkOrthogonality(problem, test, basisId, 'test');
fprintf('PASSED\n');

fprintf('\n========== LOD PML global corrector tests PASSED ==========\n');


function checkOrthogonality(problem, out, basisId, side)
[A, ~] = problem.form.global();
fineFree = problem.dof.fineFree;
C = lodMomentGlobalConstraints(problem.moment.rows, fineFree, ...
    problem.dof.coarseFree, struct('constraintTolerance', 1e-12));
Z = null(full(C'));
assert(~isempty(Z), 'Constrained global kernel is empty.');
P = problem.transfer();
v = P(fineFree, basisId) - out.correctors(fineFree, 1);
if strcmp(side, 'trial')
    residual = norm(Z' * A(fineFree, fineFree) * v) / max(1, norm(v));
else
    residual = norm(v' * A(fineFree, fineFree) * Z) / max(1, norm(v));
end
constraintResidual = norm(C' * out.correctors(fineFree, 1)) / ...
    max(1, norm(out.correctors(fineFree, 1)));
assert(residual < 1e-10, '%s variational residual %.3e too large.', side, residual);
assert(constraintResidual < 1e-10, ...
    '%s constraint residual %.3e too large.', side, constraintResidual);
end
