% VERIFY_AS_OAS_POISSON2D  AS/OAS verification and comparison for Poisson.

fprintf('========== AS/OAS Poisson 2D Verification ==========\n\n');
fprintf('Model: -Delta u = 2*pi^2*sin(pi*x)*sin(pi*y), u=0 on [0,1]^2\n');
fprintf('AS reference: kappa <= C (1 + H/delta) for overlapping one-level AS.\n');
fprintf('OAS reference: optimized Robin transmission improves Schwarz convergence.\n\n');

u_exact = @(x, y) sin(pi*x) .* sin(pi*y);
f_rhs = @(x, y) 2*pi^2*u_exact(x, y);

%% AS checkerboard condition number vs H/delta
fprintf('Table 1: One-level AS, checkerboard 4x4, exact eig condition number\n');
h = 1/24;
nGrid = [4, 4];
H = 1 / nGrid(1);
[node, elem, bd] = squaremesh([0, 1, 0, 1], h);
[A_ff, b_f, freeNodes] = poissonReducedSystem(node, elem, bd, f_rhs);

fprintf('%-8s %-8s %-8s %-10s %-10s %-8s\n', ...
    'delta/h', 'H/delta', 'Nfree', 'lambdaMin', 'lambdaMax', 'kappa');
fprintf('%s\n', repmat('-', 1, 62));

for m = [1, 2, 3, 6]
    delta = m * h;
    parts = partitionMesh2D(node, elem, bd, nGrid, 'overlap', delta);
    asPrec = additiveSchwarz(A_ff, parts, freeNodes);
    [kappa, lamMin, lamMax] = preconditionedConditionNumber(A_ff, asPrec);
    fprintf('%-8d %-8.2f %-8d %-10.3e %-10.3e %-8.2f\n', ...
        m, H/delta, length(freeNodes), lamMin, lamMax, kappa);
end

%% AS vs OAS Krylov comparison
fprintf('\nTable 2: AS(PCG) vs OAS(GMRES), checkerboard 4x4\n');
h = 1/32;
H = 1 / nGrid(1);
[node, elem, bd] = squaremesh([0, 1, 0, 1], h);
[A_ff, b_f, freeNodes] = poissonReducedSystem(node, elem, bd, f_rhs);
alphaVals = [0.5, 1, 2] * pi / H;

fprintf('%-8s %-8s %-10s %-10s %-10s %-10s\n', ...
    'delta/h', 'H/delta', 'AS_it', 'AS_relres', 'OAS_it', 'best_alpha');
fprintf('%s\n', repmat('-', 1, 70));

for m = [1, 2, 4]
    delta = m * h;
    parts = partitionMesh2D(node, elem, bd, nGrid, 'overlap', delta);

    asPrec = additiveSchwarz(A_ff, parts, freeNodes);
    [~, flagAS, relAS, iterAS] = pcg(A_ff, b_f, 1e-8, 200, asPrec);
    assert(flagAS == 0, 'AS-PCG failed for delta/h=%d', m);

    bestIt = inf;
    bestAlpha = NaN;
    for alpha = alphaVals
        oasPrec = optimizedAdditiveSchwarzPoisson2D(node, elem, bd, parts, freeNodes, alpha);
        [~, flagOAS, relOAS, iterOAS] = gmres(A_ff, b_f, [], 1e-8, 200, oasPrec);
        assert(flagOAS == 0, 'OAS-GMRES failed for delta/h=%d alpha=%.3g relres=%.3e', ...
            m, alpha, relOAS);
        itOAS = iterOAS(2);
        if itOAS < bestIt
            bestIt = itOAS;
            bestAlpha = alpha;
        end
    end

    fprintf('%-8d %-8.2f %-10d %-10.2e %-10d %-10.3g\n', ...
        m, H/delta, iterAS, relAS, bestIt, bestAlpha);
end

fprintf('\n========== AS/OAS Poisson tests PASSED ==========\n');


function [A_ff, b_f, freeNodes] = poissonReducedSystem(node, elem, bd, f_rhs)
A = assembleStiffness2D(node, elem);
M = assembleMass2D(node, elem);
b = M * f_rhs(node(:,1), node(:,2));
bdNodes = getBoundaryNodes2D(elem, bd);
freeNodes = setdiff(1:size(node, 1), bdNodes)';
A_ff = A(freeNodes, freeNodes);
b_f = b(freeNodes);
end


function [kappa, lamMin, lamMax] = preconditionedConditionNumber(A, applyPrecon)
n = size(A, 1);
MA = zeros(n, n);
for j = 1:n
    MA(:, j) = applyPrecon(A(:, j));
end

MA = (MA + MA') / 2;
eigVals = eig(MA);
eigVals = real(eigVals);
eigVals = eigVals(eigVals > 1e-10);
lamMin = min(eigVals);
lamMax = max(eigVals);
kappa = lamMax / lamMin;
end
