function fine = buildPMLLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, pml, opts)
% BUILDPMLLXZZFINESPACEHELMHOLTZ2D  P1 divergence-form PML fine algebra for LXZZ.

if nargin < 6 || isempty(opts), opts = struct(); end
opts = localOptions(opts);
if opts.degree ~= 1
    error('buildPMLLxzzFineSpaceHelmholtz2D:p1Only', ...
        'The PML LXZZ fine-space builder currently supports P1 triangles only.');
end

asmOpts = struct('quadOrder', opts.quadOrder);
[Afull, bfull, freeDof, bdDof, info] = assembleHelmholtzPMLDivergence2D( ...
    node, elem, k, pml, opts.source, opts.degree, asmOpts);
K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
energyFull = K + abs(k)^2 * M;

freeDof = freeDof(:);
bdDof = bdDof(:);
nFull = size(node, 1);
nFree = numel(freeDof);
fullToActive = zeros(nFull, 1);
fullToActive(freeDof) = (1:nFree).';

fine = struct();
fine.dim = 2;
fine.form = 'pml-divergence';
fine.degree = 1;
fine.node = node(freeDof, :);
fine.elem = zeros(0, 3);
fine.bdFlag = zeros(0, 3);
fine.baseNode = fine.node;
fine.baseElem = fine.elem;
fine.baseBdFlag = fine.bdFlag;
fine.A = Afull(freeDof, freeDof);
fine.b = bfull(freeDof);
fine.energy = energyFull(freeDof, freeDof);
fine.K = K(freeDof, freeDof);
fine.M = M(freeDof, freeDof);
fine.boundaryMass = sparse(nFree, nFree);
fine.pde = normalizeHelmholtzPDE(k);
fine.helmholtzInput = k;
fine.p1ToFine = speye(nFree);
fine.baseToFine = speye(nFree);
fine.N = nFree;
fine.pml = pml;
fine.pmlInfo = info;
fine.pmlAssemblyOptions = asmOpts;
fine.fullNode = node;
fine.fullElem = elem(:, 1:3);
fine.fullBdFlag = bdFlag;
fine.fullA = Afull;
fine.fullB = bfull;
fine.fullEnergy = energyFull;
fine.freeDof = freeDof;
fine.boundaryDof = bdDof;
fine.activeToFull = freeDof;
fine.fullToActive = fullToActive;

if opts.cacheEnergySolver
    fine.energySolve = energySolverHandle(fine.energy);
end
end


function opts = localOptions(opts)
defaults = struct();
defaults.degree = 1;
defaults.quadOrder = 4;
defaults.source = [];
defaults.cacheEnergySolver = false;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function solve = energySolverHandle(D)
try
    R = chol(D);
    solve = @(b) R \ (R' \ b);
catch
    [L, U, p, q] = lu(D, 'vector');
    solve = @(b) solveLU(L, U, p, q, b);
end
end


function x = solveLU(L, U, p, q, b)
x = zeros(size(b));
x(q, :) = U \ (L \ b(p, :));
end
