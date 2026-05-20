% VERIFY_ORAS_REPRODUCE_GGGS  ORAS experiments modelled on Gander-Gong-Graham-Spence.
%
%   Reports:
%     1. k-weighted H1 norm of E = I - B^{-1}A and E^N.
%     2. Richardson and GMRES iterations for strip/checkerboard partitions.
%     3. P1/P2/P3 comparisons at fixed points per wavelength.
%
%   The default cases are intentionally modest enough to run in this
%   repository's MATLAB verification suite.  The geometry/overlap choices
%   follow the papers.  Linear partition-of-unity weights are used by default.

fprintf('========== ORAS Gander-Gong-Graham-Spence Reproduction ==========\n\n');

fprintf('Table A: Strip power norms, p=2, domain (0,2N/3)x(0,1), overlap extension=1/6\n');
fprintf('%-5s %-5s %-7s %-8s %-10s %-10s %-10s\n', ...
    'Nsub', 'k', '1/h', 'DOF', '||E||', '||E^N||', 'GMRES');
fprintf('%s\n', repmat('-', 1, 72));
for nSub = [2, 4]
    for k = [10, 20]
        h = choosePowerMesh(k);
        [normE, normEN, gmIts, ndof] = powerExperiment('strip', nSub, k, h, 2);
        fprintf('%-5d %-5d %-7d %-8d %-10.3g %-10.3g %-10d\n', ...
            nSub, k, round(1/h), ndof, normE, normEN, gmIts);
    end
end

fprintf('\nTable B: Checkerboard power norms, p=2, unit square, overlap=H/4\n');
fprintf('%-8s %-5s %-7s %-8s %-10s %-10s %-10s\n', ...
    'grid', 'k', '1/h', 'DOF', '||E||', '||E^N||', 'GMRES');
fprintf('%s\n', repmat('-', 1, 75));
for gridN = [2, 3]
    for k = [10, 20]
        h = choosePowerMesh(k);
        [normE, normEN, gmIts, ndof] = powerExperiment('grid', gridN, k, h, 2);
        fprintf('%-8s %-5d %-7d %-8d %-10.3g %-10.3g %-10d\n', ...
            sprintf('%dx%d', gridN, gridN), k, round(1/h), ndof, normE, normEN, gmIts);
    end
end

fprintf('\nTable C: Strip iteration counts, domain (0,16/3)x(0,1), 8 strips\n');
fprintf('Resolution sweep: h = 2*pi/(q*k), q in {10,20,40,80}\n');
fprintf('%-5s %-3s %-4s %-7s %-8s %-12s %-10s\n', ...
    'k', 'p', 'q', '1/h', 'DOF', 'Richardson', 'GMRES');
fprintf('%s\n', repmat('-', 1, 65));
for k = [10, 20]
    for p = 1:3
        for q = resolutionFactors()
            h = chooseResolutionMesh(k, q);
            ndofEst = estimateDof('strip', 8, h, p);
            if ndofEst > maxIterationDof()
                fprintf('%-5d %-3d %-4d %-7d %-8d %-12s %-10s\n', ...
                    k, p, q, round(1/h), ndofEst, 'skip', 'skip');
                continue;
            end
            [richIts, gmIts, ndof] = iterationExperiment('strip', 8, k, h, p);
            fprintf('%-5d %-3d %-4d %-7d %-8d %-12s %-10s\n', ...
                k, p, q, round(1/h), ndof, iterString(richIts), iterString(gmIts));
        end
    end
end

fprintf('\nTable D: Checkerboard iteration counts, unit square, H~k^{-0.4}, overlap=H/4\n');
fprintf('Resolution sweep: h = 2*pi/(q*k), q in {10,20,40,80}\n');
fprintf('%-5s %-3s %-4s %-8s %-7s %-8s %-12s %-10s\n', ...
    'k', 'p', 'q', 'grid', '1/h', 'DOF', 'Richardson', 'GMRES');
fprintf('%s\n', repmat('-', 1, 78));
for k = [10, 20]
    gridN = max(2, round(k^0.4));
    for p = 1:3
        for q = resolutionFactors()
            h = chooseResolutionMesh(k, q);
            ndofEst = estimateDof('grid', gridN, h, p);
            if ndofEst > maxIterationDof()
                fprintf('%-5d %-3d %-4d %-8s %-7d %-8d %-12s %-10s\n', ...
                    k, p, q, sprintf('%dx%d', gridN, gridN), round(1/h), ndofEst, ...
                    'skip', 'skip');
                continue;
            end
            [richIts, gmIts, ndof] = iterationExperiment('grid', gridN, k, h, p);
            fprintf('%-5d %-3d %-4d %-8s %-7d %-8d %-12s %-10s\n', ...
                k, p, q, sprintf('%dx%d', gridN, gridN), round(1/h), ndof, ...
                iterString(richIts), iterString(gmIts));
        end
    end
end

fprintf('\n========== ORAS reproduction study complete ==========\n');


function h = choosePowerMesh(k)
% A feasible analogue of h ~ k^{-5/4}; capped to keep dense E affordable.
h = min(1/8, max(1/14, k^(-5/4)));
h = 1 / ceil(1/h);
end


function factors = resolutionFactors()
% Table 5 in GGS uses h from 2*pi/(10k) to 2*pi/(80k).
% Keep q=80 available but skip infeasible cases in iterationExperiment.
factors = [10, 20, 40, 80];
end


function h = chooseResolutionMesh(k, q)
hTarget = 2*pi/(q*k);
h = 1 / ceil(1/hTarget);
end


function nmax = maxIterationDof()
nmax = 12000;
end


function ndof = estimateDof(kind, nSub, h, degree)
switch kind
    case 'strip'
        L = 16/3;
        nx = ceil(L / h);
        ny = ceil(1 / h);
    case 'grid'
        nx = ceil(1 / h);
        ny = nx;
    otherwise
        error('Unknown kind: %s', kind);
end

switch degree
    case 1
        ndof = (nx+1) * (ny+1);
    case 2
        ndof = (2*nx+1) * (2*ny+1);
    case 3
        ndof = (3*nx+1) * (3*ny+1);
    otherwise
        error('Unsupported degree: %d', degree);
end
end


function [normE, normEN, gmIts, ndof] = powerExperiment(kind, nSub, k, h, degree)
[node, elem, bd, parts] = makeCase(kind, nSub, k, h, degree, 'power');
[A, b, D] = helmholtzCase(node, elem, bd, k, degree);
ap = orasHelmholtz(node, elem, bd, k, parts, degree);

ndof = size(A, 1);
E = errorPropagationMatrix(A, ap);
N = length(parts);
normE = weightedOperatorNorm(E, D);
EN = E^N;
normEN = weightedOperatorNorm(EN, D);

[~, flag, ~, iter] = gmres(A, b, [], 1e-6, 150, ap);
assert(flag == 0, 'GMRES failed in power experiment (%s, N=%d, k=%g).', kind, N, k);
gmIts = iter(2);
end


function [richIts, gmIts, ndof] = iterationExperiment(kind, nSub, k, h, degree)
[node, elem, bd, parts] = makeCase(kind, nSub, k, h, degree, 'iteration');
[A, b] = helmholtzCase(node, elem, bd, k, degree);
ap = orasHelmholtz(node, elem, bd, k, parts, degree);
ndof = size(A, 1);

richIts = richardsonIterations(A, b, ap, 1e-6, 120);
[~, flag, ~, iter] = gmres(A, b, [], 1e-6, 120, ap);
assert(flag == 0, 'GMRES failed in iteration experiment (%s, k=%g, p=%d).', kind, k, degree);
gmIts = iter(2);
end


function [node, elem, bd, parts] = makeCase(kind, nSub, k, h, degree, mode)
switch kind
    case 'strip'
        if strcmp(mode, 'power')
            L = 2*nSub/3;
            overlap = 1/6;
        else
            L = 16/3;
            overlap = 1/4;
            nSub = 8;
        end
        [node, elem, bd] = squaremesh([0, L, 0, 1], h);
        parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', overlap);
        parts = linearPartitionOfUnity2D(parts, [0, L, 0, 1], [nSub, 1], overlap);
    case 'grid'
        [node, elem, bd] = squaremesh([0, 1, 0, 1], h);
        if strcmp(mode, 'power')
            overlap = 1/(4*nSub);
        else
            H = 1/nSub;
            overlap = H/4;
        end
        parts = partitionMesh2D(node, elem, bd, [nSub, nSub], 'overlap', overlap);
        parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [nSub, nSub], overlap);
    otherwise
        error('Unknown ORAS case: %s', kind);
end
end


function [A, b, D] = helmholtzCase(node, elem, bd, k, degree)
if degree == 1
    nodeH = node;
    elemH = elem;
else
    [nodeH, elemH] = extendMesh2D(node, elem, degree);
end

[A, ~] = assembleHelmholtz2D(nodeH, elemH, bd, k, 0, 0, degree);
uEx = exp(1i*k*nodeH(:,1));
b = A * uEx;

if nargout > 2
    K = assembleStiffness2D(nodeH, elemH, degree);
    M = assembleMass2D(nodeH, elemH, degree);
    D = K + k^2*M;
    D = (D + D') / 2;
end
end


function E = errorPropagationMatrix(A, ap)
n = size(A, 1);
E = eye(n);
for j = 1:n
    E(:, j) = E(:, j) - ap(A(:, j));
end
end


function opNorm = weightedOperatorNorm(E, D)
S = E' * D * E;
S = (S + S') / 2;
ev = eig(full(S), full(D));
opNorm = sqrt(max(real(ev)));
end


function its = richardsonIterations(A, b, ap, tol, maxIter)
u = zeros(size(b));
r0 = norm(b);
its = maxIter + 1;
for it = 1:maxIter
    r = b - A*u;
    if norm(r) / r0 < tol
        its = it - 1;
        return;
    end
    u = u + ap(r);
end
if norm(b - A*u) / r0 < tol
    its = maxIter;
end
end


function s = iterString(its)
if isinf(its)
    s = 'skip';
elseif its > 120
    s = '>120';
else
    s = sprintf('%d', its);
end
end
