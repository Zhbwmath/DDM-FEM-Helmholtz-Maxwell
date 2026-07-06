function lod = buildLODHelmholtzPML2D(nodeH, elemH, bdH, nodeh, elemh, bdh, k, pml, f, opts)
% BUILDLODHELMHOLTZPML2D  Build L2-moment LOD for 2D Helmholtz PML.

if nargin < 9, f = []; end
if nargin < 10 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'solveCoarse') || isempty(opts.solveCoarse), opts.solveCoarse = true; end
if ~isfield(opts, 'correctorSide') || isempty(opts.correctorSide), opts.correctorSide = 'both'; end

problem = helmholtzPMLLODProblem2D(nodeH, elemH, bdH, nodeh, elemh, bdh, ...
    k, pml, f, opts);
lodOpts = opts;
lodOpts.solveCoarse = false;
lod = buildLOD(nodeH, elemH, nodeh, elemh, problem, lodOpts);

coarseFree = problem.dof.coarseFree(:);
fineFree = problem.dof.fineFree(:);
P = lod.basis.trial(:, coarseFree);
Pte = lod.basis.test(:, coarseFree);
A = lod.system.A;
b = lod.system.b;
AH = Pte' * A * P;
bH = Pte' * b;

lod.system.AH = AH;
lod.system.bH = bH;
lod.system.coarseFree = coarseFree;
lod.problem = struct('name', 'HelmholtzPML2D', 'k', k, ...
    'form', 'divergence PML with L2 moment constraints', 'pml', pml);
lod.moment = problem.moment;
lod.dof = problem.dof;
lod.solution = struct('xH', [], 'uh', [], 'fine', [], ...
    'relativePhysicalEnergyError', NaN, 'relativePhysicalL2Error', NaN);

if opts.solveCoarse
    xH = AH \ bH;
    uh = P * xH;
    uf = zeros(size(nodeh, 1), 1);
    uf(fineFree) = A(fineFree, fineFree) \ b(fineFree);
    [Ephys, Mphys] = physicalEnergyMatrices(nodeh, elemh, k, pml);
    e = uf - uh;
    lod.solution.xH = xH;
    lod.solution.uh = uh;
    lod.solution.fine = uf;
    lod.solution.relativePhysicalEnergyError = energyNorm(e, Ephys) / ...
        max(1, energyNorm(uf, Ephys));
    lod.solution.relativePhysicalL2Error = energyNorm(e, Mphys) / ...
        max(1, energyNorm(uf, Mphys));
end
end


function [E, M] = physicalEnergyMatrices(node, elem, k, pml)
physElem = physicalElements(node, elem, pml.physicalBox);
K = assembleStiffness2D(node, elem(physElem, :), 1);
M = assembleMass2D(node, elem(physElem, :), 1);
E = K + k^2 * M;
end


function isPhys = physicalElements(node, elem, box)
centroid = (node(elem(:,1), :) + node(elem(:,2), :) + node(elem(:,3), :)) / 3;
tol = 100 * eps(max(1, max(abs(box))));
isPhys = centroid(:,1) >= box(1) - tol & centroid(:,1) <= box(2) + tol & ...
    centroid(:,2) >= box(3) - tol & centroid(:,2) <= box(4) + tol;
end


function nrm = energyNorm(v, E)
val = real(v' * E * v);
nrm = sqrt(max(0, val));
end
