% VERIFY_LXZZ25_ADJOINT_GMRES_VARIANTS  Compare reference adjoint and Dk-GMRES.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

tol = 1e-6;
maxit = envNumber('LXZZ25_VARIANT_MAXIT', 40);
kValues = parseNumberList(envString('LXZZ25_VARIANT_KVALUES', '16,32'));
variants = {'dirichlet', 'impedance'};
variantFilter = envString('LXZZ25_VARIANT_TYPES', '');
if ~isempty(variantFilter)
    variants = cellstr(split(string(variantFilter), ',')).';
end
methodFilter = lower(envString('LXZZ25_VARIANT_METHODS', 'reference,dk'));

fprintf('========== LXZZ25 adjoint / GMRES variant check ==========\n\n');
fprintf('Rows use exact Section 5.1/5.2 small-case parameters.\n');
fprintf('euclidean-reference: MATLAB GMRES with reference Q0^* correction.\n');
fprintf('Dk-model: unrestarted Arnoldi/least-squares in D_k inner product with D_k adjoint.\n\n');

for ik = 1:numel(kValues)
    k = kValues(ik);
    for iv = 1:numel(variants)
        variant = variants{iv};
        c = smallCase(k, variant);
        fprintf('k=%g %-10s h=1/%d H=1/%d Hsub=1/%d m=%d\n', ...
            k, variant, c.hInv, c.HInv, c.HsubInv, c.m);

        [node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1 / c.hInv);
        [nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1 / c.HInv);
        parts = coarseHatPartition2D(node, elem, bdFlag, c.Hsub / 2);
        b = assemblePlaneWaveBoundaryLoadP1(node, elem, bdFlag, k);

        lodOpts = struct('oversampling', c.m, 'solveCoarse', false, ...
            'solverMode', 'direct', 'useParfor', false);

        if contains(methodFilter, 'reference')
            optsRef = struct('variant', variant, 'coarseType', 'lod', ...
                'lodOptions', lodOpts, 'solverMode', 'lu', 'useParfor', false, ...
                'adjointType', 'reference');
            preRef = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, ...
                parts, nodeH, elemH, bdH, optsRef);
            [~, flagRef, relRef, iterRef, resvecRef] = gmres(preRef.A, b, [], tol, maxit, @preRef.applyResidual);
            fprintf('  euclidean-reference: flag=%d iter=%d relres=%.3e first=%s\n', ...
                flagRef, gmresIterationCount(iterRef), relRef, shortResidualHistory(resvecRef));
        end

        if contains(methodFilter, 'dk')
            optsDk = struct('variant', variant, 'coarseType', 'lod', ...
                'lodOptions', lodOpts, 'solverMode', 'lu', 'useParfor', false, ...
                'adjointType', 'energy');
            preDk = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, ...
                parts, nodeH, elemH, bdH, optsDk);
            outDk = runWeightedModelGmres(preDk, b, tol, maxit);
            fprintf('  Dk-model:            flag=%d iter=%d relres=%.3e precondDrel=%.3e first=%s\n', ...
                outDk.flag, outDk.iter, outDk.originalRelres, outDk.precondDRelres, ...
                shortNumericHistory(outDk.originalResidualHistory));
            fprintf('    timing: Aop %.3fs, D %.3fs, orth %.3fs, ls %.3fs, monitor %.3fs\n', ...
                outDk.timing.operatorApply, outDk.timing.weightApply, ...
                outDk.timing.orthogonalization, outDk.timing.leastSquares, outDk.timing.monitor);
        end
        fprintf('\n');
    end
end

fprintf('========== variant check complete ==========\n');


function c = smallCase(k, variant)
if strcmpi(variant, 'dirichlet')
    targetHInv = k;
    targetHsubInv = k / 2;
else
    targetHInv = k;
    targetHsubInv = k / 4;
end
raw = ceil(k^(3/2));
hInv = alignFineInv(raw, [targetHInv, targetHsubInv]);
c = struct();
c.hInv = hInv;
c.HInv = compatibleDivisor(hInv, targetHInv);
c.HsubInv = compatibleDivisor(hInv, targetHsubInv);
c.Hsub = 1 / c.HsubInv;
c.m = round(log2(k)) - 1;
end


function out = runWeightedModelGmres(pre, b, tol, maxit)
bnorm = norm(b);
if bnorm == 0, bnorm = 1; end
rhs = pre.applyResidual(b);
opts = struct();
opts.monitor = @(x,~,~) norm(b - pre.A * x) / bnorm;
opts.collectTiming = true;
[~, flag, relres, iter, ~, info] = weightedGmres(@pre.apply, rhs, pre.energy, tol, maxit, opts);
out = struct();
out.flag = flag;
out.iter = iter;
out.originalRelres = info.monitorValues(end);
out.precondDRelres = relres;
out.originalResidualHistory = info.monitorValues;
out.timing = info.timing;
end


function b = assemblePlaneWaveBoundaryLoadP1(node, elem, bdFlag, k)
N = size(node, 1);
edgeVertex = [2 3; 3 1; 1 2];
b = zeros(N, 1);
for t = 1:size(elem, 1)
    for e = 1:3
        if bdFlag(t, e) ~= 1, continue; end
        va = elem(t, edgeVertex(e, 1));
        vb = elem(t, edgeVertex(e, 2));
        mid = 0.5 * (node(va, :) + node(vb, :));
        normal = squareBoundaryNormal(mid);
        ga = planeWaveBoundaryValue(node(va, :), normal, k);
        gb = planeWaveBoundaryValue(node(vb, :), normal, k);
        L = norm(node(vb, :) - node(va, :));
        b([va; vb]) = b([va; vb]) + L / 6 * [2, 1; 1, 2] * [ga; gb];
    end
end
end


function n = squareBoundaryNormal(x)
tol = 1e-12;
if abs(x(1) - 1) < tol
    n = [1, 0];
elseif abs(x(2)) < tol
    n = [0, -1];
elseif abs(x(1)) < tol
    n = [-1, 0];
elseif abs(x(2) - 1) < tol
    n = [0, 1];
else
    error('verify_lxzz25:normal', 'Boundary midpoint was not on the unit square.');
end
end


function g = planeWaveBoundaryValue(x, n, k)
d = [1 / sqrt(2), 1 / sqrt(2)];
u = exp(1i * k * (x(:, 1) * d(1) + x(:, 2) * d(2)));
g = 1i * k * (d * n.' - 1) .* u;
end


function n = gmresIterationCount(iter)
if isempty(iter)
    n = NaN;
elseif numel(iter) == 1
    n = iter;
else
    n = iter(2);
end
end


function s = shortResidualHistory(resvec)
rr = resvec ./ resvec(1);
s = shortNumericHistory(rr);
end


function s = shortNumericHistory(v)
n = min(numel(v), 10);
parts = cell(1, n);
for i = 1:n
    parts{i} = sprintf('%d:%.1e', i - 1, v(i));
end
s = strjoin(parts, ' ');
end


function n = alignFineInv(raw, targetInvs)
targetInvs = round(targetInvs(isfinite(targetInvs) & targetInvs > 0));
base = 1;
for i = 1:numel(targetInvs)
    base = lcm(base, max(1, targetInvs(i)));
end
n = ceil(raw / base) * base;
end


function d = compatibleDivisor(n, target)
target = max(1, target);
divs = divisorsInt(n);
[~, idx] = min(abs(divs - target));
d = divs(idx);
end


function d = divisorsInt(n)
d = [];
for q = 1:floor(sqrt(n))
    if mod(n, q) == 0
        d(end+1) = q; %#ok<AGROW>
        if q ~= n / q
            d(end+1) = n / q; %#ok<AGROW>
        end
    end
end
d = sort(d);
end


function v = envNumber(name, defaultValue)
s = getenv(name);
if isempty(s)
    v = defaultValue;
else
    v = str2double(s);
    if isnan(v), v = defaultValue; end
end
end


function v = envString(name, defaultValue)
s = getenv(name);
if isempty(s)
    v = defaultValue;
else
    v = s;
end
end


function values = parseNumberList(s)
parts = split(string(s), ',');
values = zeros(1, numel(parts));
for i = 1:numel(parts)
    values(i) = str2double(parts(i));
end
values = values(isfinite(values));
end
