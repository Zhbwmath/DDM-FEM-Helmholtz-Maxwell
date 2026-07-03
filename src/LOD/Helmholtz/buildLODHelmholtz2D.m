function lod = buildLODHelmholtz2D(nodeH, elemH, bdH, nodeh, elemh, bdh, k, f, g, opts)
% BUILDLODHELMHOLTZ2D  Build Petrov-Galerkin LOD for 2D Helmholtz.

if nargin < 8, f = []; end
if nargin < 9, g = []; end
if nargin < 10 || isempty(opts), opts = struct(); end

kInput = helmholtzLODInput(k);
problem = helmholtzLODProblem2D(nodeH, elemH, bdH, nodeh, elemh, bdh, kInput, f, g);
lod = buildLOD(nodeH, elemH, nodeh, elemh, problem, opts);
lod.problem = struct('name', 'Helmholtz2D', 'k', k, ...
    'form', 'grad-grad - k^2 mass - i k boundary mass');

K = assembleStiffness2D(nodeh, elemh, 1);
if isnumeric(kInput) && isscalar(kInput)
    M = assembleMass2D(nodeh, elemh, 1);
    lod.system.energy = K + (kInput^2) * M;
else
    pde = normalizeHelmholtzPDE(kInput);
    qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
    lod.system.energy = K + assembleWeightedMass2D(nodeh, elemh, 1, qfun);
end

if ~isempty(lod.solution.uh)
    uf = lod.system.A \ lod.system.b;
    e = uf - lod.solution.uh;
    lod.solution.fine = uf;
    lod.solution.relativeEnergyError = sqrt(real(e' * lod.system.energy * e)) / ...
        max(1, sqrt(real(uf' * lod.system.energy * uf)));
end
end


function input = helmholtzLODInput(k)
if isnumeric(k) && isscalar(k)
    input = k;
else
    input = normalizeHelmholtzPDE(k);
end
end
