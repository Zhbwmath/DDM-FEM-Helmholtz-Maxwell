function lod = buildLODHelmholtz3D(nodeH, elemH, bdH, nodeh, elemh, bdh, k, f, g, opts)
% BUILDLODHELMHOLTZ3D  Build Petrov-Galerkin LOD for 3D Helmholtz.

if nargin < 8, f = []; end
if nargin < 9, g = []; end
if nargin < 10 || isempty(opts), opts = struct(); end

problem = helmholtzLODProblem3D(nodeH, elemH, bdH, nodeh, elemh, bdh, k, f, g);
lod = buildLOD(nodeH, elemH, nodeh, elemh, problem, opts);
lod.problem = struct('name', 'Helmholtz3D', 'k', k, ...
    'form', 'grad-grad - k^2 mass - i k boundary mass');

K = assembleStiffness3D(nodeh, elemh, 1);
M = assembleMass3D(nodeh, elemh, 1);
lod.system.energy = K + (k^2) * M;

if ~isempty(lod.solution.uh)
    uf = lod.system.A \ lod.system.b;
    e = uf - lod.solution.uh;
    lod.solution.fine = uf;
    lod.solution.relativeEnergyError = sqrt(real(e' * lod.system.energy * e)) / ...
        max(1, sqrt(real(uf' * lod.system.energy * uf)));
end
end
