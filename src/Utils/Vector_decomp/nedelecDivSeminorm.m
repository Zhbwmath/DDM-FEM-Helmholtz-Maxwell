function [seminorm, info] = nedelecDivSeminorm(node, elem, u, opts)
% NEDELECDIVSEMINORM  Compute the weak P1-gradient discrete divergence seminorm.

if nargin < 4 || isempty(opts), opts = struct(); end

[G, edge] = nedelecGradientMatrix(node, elem, opts);
NE = size(G, 1);
N = size(G, 2);
dim = size(node, 2);

if isvector(u) && numel(u) == NE
    u = u(:);
end
if size(u, 1) ~= NE
    error('nedelecDivSeminorm:badVectorSize', ...
        'The Nedelec coefficient array must have %d rows.', NE);
end

if isfield(opts, 'massMatrix') && ~isempty(opts.massMatrix)
    M = opts.massMatrix;
else
    switch dim
        case 2
            M = assembleNedMass2D(node, elem);
        case 3
            M = assembleNedMass3D(node, elem);
        otherwise
            error('nedelecDivSeminorm:badDimension', ...
                'Only 2D triangular and 3D tetrahedral meshes are supported.');
    end
end
if ~isequal(size(M), [NE, NE])
    error('nedelecDivSeminorm:badMassSize', ...
        'The Nedelec mass matrix must be %d-by-%d.', NE, NE);
end

if isfield(opts, 'edgeDofs') && ~isempty(opts.edgeDofs)
    edgeDofs = opts.edgeDofs(:);
else
    edgeDofs = (1:NE).';
end
if isfield(opts, 'scalarDofs') && ~isempty(opts.scalarDofs)
    scalarDofs = opts.scalarDofs(:);
else
    scalarDofs = (1:N).';
end
if isfield(opts, 'fixedScalarDofs')
    fixedScalarDofs = opts.fixedScalarDofs(:);
else
    fixedScalarDofs = [];
end
if isempty(fixedScalarDofs) && numel(scalarDofs) == N
    fixedScalarDofs = scalarDofs(1);
end

solveDofs = setdiff(scalarDofs, fixedScalarDofs, 'stable');
Mloc = M(edgeDofs, edgeDofs);
Gsolve = G(edgeDofs, solveDofs);
r = Gsolve' * (Mloc * u(edgeDofs, :));

if isempty(solveDofs)
    y = zeros(0, size(u, 2), class(u));
    sq = zeros(1, size(u, 2), class(u));
else
    H = Gsolve' * (Mloc * Gsolve);
    y = H \ r;
    sq = real(sum(conj(r) .* y, 1));
end
seminorm = sqrt(max(sq, 0));

info = struct();
info.gradientMatrix = G;
info.massMatrix = M;
info.edge = edge;
info.edgeDofs = edgeDofs;
info.scalarDofs = scalarDofs;
info.fixedScalarDofs = fixedScalarDofs;
info.solveScalarDofs = solveDofs;
info.weakDivergenceResidual = r;
info.weakDivergenceResidualNorm = norm(r, 'fro');
info.dualPotential = y;
info.squaredSeminorm = sq;
info.relativeResidual = info.weakDivergenceResidualNorm / ...
    max(norm(Gsolve, 'fro') * norm(Mloc * u(edgeDofs, :), 'fro'), eps);
end
