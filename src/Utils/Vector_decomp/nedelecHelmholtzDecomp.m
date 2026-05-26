function [uGrad, uComp, phi, info] = nedelecHelmholtzDecomp(node, elem, u, opts)
% NEDELECHELMHOLTZDECOMP  Split u_h = grad phi_h + z_h with G'*M*z_h = 0.

if nargin < 4 || isempty(opts), opts = struct(); end

[G, edge] = nedelecGradientMatrix(node, elem, opts);
NE = size(G, 1);
N = size(G, 2);
dim = size(node, 2);

if isvector(u) && numel(u) == NE
    u = u(:);
end
if size(u, 1) ~= NE
    error('nedelecHelmholtzDecomp:badVectorSize', ...
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
    end
end
if ~isequal(size(M), [NE, NE])
    error('nedelecHelmholtzDecomp:badMassSize', ...
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
phi = zeros(N, size(u, 2), class(u));

Mloc = M(edgeDofs, edgeDofs);
Gsolve = G(edgeDofs, solveDofs);
uloc = u(edgeDofs, :);

if ~isempty(solveDofs)
    H = Gsolve' * (Mloc * Gsolve);
    rhs = Gsolve' * (Mloc * uloc);
    alpha = H \ rhs;
    phi(solveDofs, :) = alpha;
    normalResidual = norm(H * alpha - rhs, 'fro') / max(norm(rhs, 'fro'), eps);
else
    H = sparse(0, 0);
    rhs = zeros(0, size(u, 2), class(u));
    normalResidual = 0;
end

uGrad = G * phi;
uComp = u - uGrad;

orth = Gsolve' * (Mloc * uComp(edgeDofs, :));
info = struct();
info.gradientMatrix = G;
info.massMatrix = M;
info.normalMatrix = H;
info.normalRhs = rhs;
info.edge = edge;
info.edgeDofs = edgeDofs;
info.scalarDofs = scalarDofs;
info.fixedScalarDofs = fixedScalarDofs;
info.solveScalarDofs = solveDofs;
info.orthogonalityNorm = norm(orth, 'fro');
info.relativeOrthogonality = info.orthogonalityNorm / ...
    max(norm(Gsolve, 'fro') * norm(Mloc * uComp(edgeDofs, :), 'fro'), eps);
info.normalResidual = normalResidual;
info.reconstructionError = norm(u - uGrad - uComp, 'fro') / max(norm(u, 'fro'), eps);
end
