function [q, lambda, info] = lodSolveConstrainedSaddle(A, C, R, opts)
% LODSOLVECONSTRAINEDSADDLE  Solve [A C; C' 0][q;lambda]=[R;0].

if nargin < 4 || isempty(opts), opts = struct(); end
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
    lambda = zeros(0, size(R, 2));
    relRes = norm(A * q - R, 'fro') / max(1, norm(R, 'fro'));
    conRes = 0;
else
    saddle = [A, C; C', sparse(nc, nc)];
    rhs = [R; sparse(nc, size(R, 2))];
    switch lower(opts.solverMode)
        case 'direct'
            sol = saddle \ rhs;
        case {'lu', 'factorized'}
            D = decomposition(saddle, 'lu');
            sol = D \ rhs;
        otherwise
            error('lodSolveConstrainedSaddle:solverMode', ...
                'Unknown solverMode "%s".', opts.solverMode);
    end
    q = sol(1:n, :);
    lambda = sol(n + (1:nc), :);
    relRes = norm(saddle * sol - rhs, 'fro') / max(1, norm(rhs, 'fro'));
    conRes = norm(C' * q, 'fro') / max(1, norm(q, 'fro'));
end

info = struct('relativeResidual', relRes, ...
    'constraintResidual', conRes, ...
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
