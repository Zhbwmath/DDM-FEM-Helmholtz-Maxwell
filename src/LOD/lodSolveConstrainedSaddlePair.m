function [q, qStar, lambda, lambdaStar, info, infoStar] = ...
    lodSolveConstrainedSaddlePair(A, C, R, Rstar, opts)
% LODSOLVECONSTRAINEDSADDLEPAIR  Solve primal/adjoint constrained systems.
%
% Reuses one factorization of [A C; C' 0]. The adjoint constrained matrix is
% the conjugate transpose of this saddle matrix.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'constraintTolerance') || isempty(opts.constraintTolerance)
    opts.constraintTolerance = 1e-12;
end
if ~isfield(opts, 'dropDependentConstraints') || isempty(opts.dropDependentConstraints)
    opts.dropDependentConstraints = true;
end
if ~isfield(opts, 'solverMode') || isempty(opts.solverMode)
    opts.solverMode = 'direct';
end

n = size(A, 1);
if isempty(C)
    C = sparse(n, 0);
end

keptColumns = 1:size(C, 2);
if opts.dropDependentConstraints && ~isempty(C)
    [C, keptColumns] = independentConstraintColumns(C, opts.constraintTolerance);
end

nc = size(C, 2);
if nc == 0
    q = A \ R;
    qStar = A' \ Rstar;
    lambda = zeros(0, size(R, 2));
    lambdaStar = zeros(0, size(Rstar, 2));
    relRes = norm(A * q - R, 'fro') / max(1, norm(R, 'fro'));
    relResStar = norm(A' * qStar - Rstar, 'fro') / max(1, norm(Rstar, 'fro'));
    conRes = 0;
    conResStar = 0;
else
    saddle = [A, C; C', sparse(nc, nc)];
    rhs = [R; sparse(nc, size(R, 2))];
    rhsStar = [Rstar; sparse(nc, size(Rstar, 2))];
    switch lower(opts.solverMode)
        case 'direct'
            sol = saddle \ rhs;
            solStar = saddle' \ rhsStar;
        case {'lu', 'factorized'}
            D = decomposition(saddle, 'lu');
            sol = D \ rhs;
            solStar = D' \ rhsStar;
        otherwise
            error('lodSolveConstrainedSaddlePair:solverMode', ...
                'Unknown solverMode "%s".', opts.solverMode);
    end
    q = sol(1:n, :);
    lambda = sol(n + (1:nc), :);
    qStar = solStar(1:n, :);
    lambdaStar = solStar(n + (1:nc), :);
    relRes = norm(saddle * sol - rhs, 'fro') / max(1, norm(rhs, 'fro'));
    relResStar = norm(saddle' * solStar - rhsStar, 'fro') / max(1, norm(rhsStar, 'fro'));
    conRes = norm(C' * q, 'fro') / max(1, norm(q, 'fro'));
    conResStar = norm(C' * qStar, 'fro') / max(1, norm(qStar, 'fro'));
end

info = struct('relativeResidual', relRes, ...
    'constraintResidual', conRes, ...
    'keptConstraintColumns', keptColumns, ...
    'numConstraints', nc);
infoStar = struct('relativeResidual', relResStar, ...
    'constraintResidual', conResStar, ...
    'keptConstraintColumns', keptColumns, ...
    'numConstraints', nc);
end


function [Ckeep, kept] = independentConstraintColumns(C, tol)
colNorm = sqrt(sum(abs(C).^2, 1));
nonzero = find(colNorm > tol * max(1, max(colNorm)));
if isempty(nonzero)
    Ckeep = sparse(size(C, 1), 0);
    kept = [];
    return;
end

Cdense = full(C(:, nonzero));
[~, R, E] = qr(Cdense, 0);
d = abs(diag(R));
rankC = nnz(d > tol * max(1, d(1)));
if rankC == 0
    Ckeep = sparse(size(C, 1), 0);
    kept = [];
    return;
end
kept = sort(nonzero(E(1:rankC)));
Ckeep = C(:, kept);
end
