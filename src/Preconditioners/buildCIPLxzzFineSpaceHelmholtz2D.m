function fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, opts)
% BUILDCIPLXZZFINESPACEHELMHOLTZ2D  CIP fine-space data for LXZZ hybrids.

if nargin < 5 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

baseElem = elem(:, 1:3);
baseNode = node(1:max(baseElem(:)), :);
if opts.degree == 1
    fineNode = baseNode;
    fineElem = baseElem;
    p1ToFine = speye(size(baseNode, 1));
elseif opts.degree == 2
    if size(elem, 2) == 3
        [fineNode, fineElem] = extendMesh2D(baseNode, baseElem, 2);
    else
        fineNode = node;
        fineElem = elem;
    end
    p1ToFine = prolongate_P1_P2(baseNode, baseElem);
elseif opts.degree == 3
    if size(elem, 2) == 3
        [fineNode, fineElem] = extendMesh2D(baseNode, baseElem, 3);
    else
        fineNode = node;
        fineElem = elem;
    end
    p1ToFine = prolongate_P1_P3(baseNode, baseElem);
else
    error('buildCIPLxzzFineSpaceHelmholtz2D:degree', ...
        'Only degree 1, 2, and 3 are supported.');
end

helmholtzInput = helmholtzMatrixInput(k);
pde = normalizeHelmholtzPDE(helmholtzInput);
K = assembleStiffness2D(fineNode, fineElem, opts.degree);
M = assembleMass2D(fineNode, fineElem, opts.degree);
Mb = assembleBoundaryMass2D(fineNode, fineElem, bdFlag, opts.degree);
[A, ~, C] = assembleHelmholtzCIP2D(fineNode, fineElem, bdFlag, ...
    helmholtzInput, [], [], opts.degree, opts.gamma, opts.cipOptions);

if isnumeric(helmholtzInput) && isscalar(helmholtzInput)
    energy = K + (abs(helmholtzInput)^2) * M;
else
    qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
    energy = K + assembleWeightedMass2D(fineNode, fineElem, ...
        opts.degree, qfun);
end

fine = struct();
fine.dim = 2;
fine.form = 'cip';
fine.degree = opts.degree;
fine.node = fineNode;
fine.elem = fineElem;
fine.bdFlag = bdFlag;
fine.baseNode = baseNode;
fine.baseElem = baseElem;
fine.baseBdFlag = bdFlag;
fine.K = K;
fine.M = M;
fine.boundaryMass = Mb;
fine.A = A;
fine.energy = energy;
fine.pde = pde;
fine.helmholtzInput = helmholtzInput;
fine.p1ToFine = p1ToFine;
fine.baseToFine = p1ToFine;
fine.N = size(fineNode, 1);
fine.cip = struct('matrix', C, 'gamma', opts.gamma, ...
    'options', opts.cipOptions, 'description', ...
    'CIP jump-stabilized Helmholtz fine form');

if opts.cacheEnergySolver
    fine.energySolve = energySolverHandle(energy);
end
end


function opts = localOptions(opts)
defaults = struct();
defaults.degree = 1;
defaults.gamma = [];
defaults.cipOptions = struct();
defaults.cacheEnergySolver = false;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function input = helmholtzMatrixInput(k)
if isnumeric(k) && isscalar(k)
    input = k;
else
    input = normalizeHelmholtzPDE(k);
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
