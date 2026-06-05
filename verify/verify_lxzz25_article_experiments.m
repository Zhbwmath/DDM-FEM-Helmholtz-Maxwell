% VERIFY_LXZZ25_ARTICLE_EXPERIMENTS  Run LXZZ25 Section 5 cases under HPC gates.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

cfg = lxzz25Config();
outDir = fullfile(repoRoot, 'verify', 'lxzz25_hybrid_lod');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cases = lxzz25Cases(cfg);
results = repmat(emptyResult(), numel(cases), 1);

fprintf('========== LXZZ25 Article Experiments ==========\n');
fprintf('Cases: %d, memory limit: %.1f GB, interactive DOF cap: %d\n\n', ...
    numel(cases), cfg.memoryLimitGB, cfg.maxRunDof);

if cfg.useParfor && cfg.runEnabled
    ensureParallelPool();
end

for i = 1:numel(cases)
    c = finalizeCase(cases(i), cfg);
    est = estimateCase(c, cfg);
    r = emptyResult();
    r = copyCaseFields(r, c, est);

    fprintf('[%02d/%02d] Table %s %s/%s k=%g %s: N=%d, estimate %.2f GB ... ', ...
        i, numel(cases), c.tableId, c.variant, c.coarseType, c.kappa, c.parameter, est.N, est.totalGB);

    if est.totalGB > cfg.memoryLimitGB
        r.status = 'blocked_memory_gt_limit';
        r.notes = sprintf('Estimate %.2f GB exceeds %.2f GB permission threshold.', ...
            est.totalGB, cfg.memoryLimitGB);
        fprintf('BLOCKED\n');
    elseif ~cfg.runEnabled
        r.status = 'estimated_only';
        r.notes = 'Execution disabled by LXZZ25_RUN=0.';
        fprintf('ESTIMATED\n');
    elseif ~cfg.runAllPermitted && (est.N > cfg.maxRunDof || est.nCoarseElem > cfg.maxRunCoarseElem)
        r.status = 'queued_runtime_cap';
        r.notes = sprintf('Below memory gate but above interactive cap: N=%d, coarse elements=%d.', ...
            est.N, est.nCoarseElem);
        fprintf('QUEUED\n');
    else
        try
            runOut = runOneCase(c, cfg);
            r = mergeRunOutput(r, runOut);
            if r.flag == 0
                r.status = 'ran';
            else
                r.status = 'ran_not_converged';
            end
            r.notes = '';
            fprintf('RAN (%d it, relres %.3e)\n', r.gmresIter, r.finalRelres);
        catch ME
            r.status = 'failed';
            r.notes = sprintf('%s: %s', ME.identifier, ME.message);
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    results(i) = r;
    writeOutputs(outDir, results(1:i), cfg);
end

writeOutputs(outDir, results, cfg);
fprintf('\n========== LXZZ25 article experiment pass complete ==========\n');
fprintf('CSV: %s\n', fullfile(outDir, 'lxzz25_article_results.csv'));
fprintf('Markdown: %s\n', fullfile(outDir, 'lxzz25_article_results.md'));


function cfg = lxzz25Config()
cfg.memoryLimitGB = envNumber('LXZZ25_MEMORY_LIMIT_GB', 300);
cfg.maxRunDof = envNumber('LXZZ25_MAX_RUN_DOF', 50000);
cfg.maxRunCoarseElem = envNumber('LXZZ25_MAX_RUN_COARSE_ELEM', 3000);
cfg.runAllPermitted = logical(envNumber('LXZZ25_RUN_ALL_PERMITTED', 0));
cfg.runEnabled = logical(envNumber('LXZZ25_RUN', 1));
cfg.useParfor = logical(envNumber('LXZZ25_PARFOR', 1));
cfg.solverMode = envString('LXZZ25_SOLVER_MODE', 'lu');
cfg.lodSolverMode = envString('LXZZ25_LOD_SOLVER_MODE', 'direct');
cfg.tol = envNumber('LXZZ25_TOL', 1e-6);
cfg.maxit = envNumber('LXZZ25_MAXIT', 100);
cfg.restart = [];
cfg.Ch = envNumber('LXZZ25_CH', 1);
cfg.CH = envNumber('LXZZ25_CH_COARSE', 1);
cfg.luFillConstant = envNumber('LXZZ25_LU_FILL', 40);
cfg.gmresBasisLength = envNumber('LXZZ25_GMRES_BASIS', 103);
cfg.timeLimitS = envNumber('LXZZ25_TIME_LIMIT_S', 7200);
cfg.adjointType = envString('LXZZ25_ADJOINT_TYPE', 'energy');
end


function cases = lxzz25Cases(cfg)
cases = repmat(baseCase(), 0, 1);
ksWave = [16, 32, 64, 128, 256, 500];
paper51 = [9, 8, 8, 8, 8, 8];
paper52 = [7, 7, 7, 7, 7, 7];
for i = 1:numel(ksWave)
    k = ksWave(i);
    m = round(log2(k)) - 1;
    cases(end+1) = makeCase('5.1', 'dirichlet', 'lod', k, 'k-derived', NaN, k, k/2, 1/k, m, ...
        sprintf('k=%g', k), paper51(i)); %#ok<AGROW>
    cases(end+1) = makeCase('5.2', 'impedance', 'lod', k, 'k-derived', NaN, k, k/4, 2/k, m, ...
        sprintf('k=%g', k), paper52(i)); %#ok<AGROW>
end

hsPow = 10:13;
paper53 = [10, 9, 9, 9];
paper54 = [9, 8, 8, 8];
for i = 1:numel(hsPow)
    hinv = 2^hsPow(i);
    cases(end+1) = makeCase('5.3', 'dirichlet', 'lod', 80, 'fixed', hinv, 80, 40, 1/80, 2, ...
        sprintf('h=2^-%d', hsPow(i)), paper53(i)); %#ok<AGROW>
    cases(end+1) = makeCase('5.4', 'impedance', 'lod', 80, 'fixed', hinv, 80, 20, 2/80, 2, ...
        sprintf('h=2^-%d', hsPow(i)), paper54(i)); %#ok<AGROW>
end

ms = [6, 5, 4, 3, 2, 1];
paper55 = [8, 8, 8, 8, 9, 11];
paper56 = [7, 7, 7, 7, 8, 10];
for i = 1:numel(ms)
    cases(end+1) = makeCase('5.5', 'dirichlet', 'lod', 128, 'k-derived', NaN, 128, 64, 1/128, ms(i), ...
        sprintf('m=%d', ms(i)), paper55(i)); %#ok<AGROW>
    cases(end+1) = makeCase('5.6', 'impedance', 'lod', 128, 'k-derived', NaN, 128, 32, 2/128, ms(i), ...
        sprintf('m=%d', ms(i)), paper56(i)); %#ok<AGROW>
end

ksMid = [40, 80, 120, 160];
paper57Lod = [9, 8; 9, 8; 9, 8; 9, 8];
paper57P1 = {'25', '23'; '55', '47'; '>100', '85'; '>100', '>100'};
for i = 1:numel(ksMid)
    k = ksMid(i);
    for j = 1:2
        if j == 1
            HsubInv = k/2; Hinv = k; label = 'Hsub=2H0';
        else
            HsubInv = k; Hinv = 2*k; label = 'Hsub=H0';
        end
        cases(end+1) = makeCase('5.7', 'dirichlet', 'lod', k, 'k-derived', NaN, Hinv, HsubInv, 1/Hinv, 2, ...
            label, paper57Lod(i, j)); %#ok<AGROW>
        cases(end+1) = makeCase('5.7', 'dirichlet', 'p1', k, 'k-derived', NaN, Hinv, HsubInv, 1/Hinv, 2, ...
            label, paper57P1{i, j}); %#ok<AGROW>
    end
end

paper58Lod = [7, 6; 7, 6; 7, 6; 6, 6];
paper58P1 = {'26', '21'; '49', '43'; '89', '77'; '>100', '>100'};
for i = 1:numel(ksMid)
    k = ksMid(i);
    for j = 1:2
        if j == 1
            HsubInv = k/4; delta = 2/k; label = 'delta=2H0';
        else
            HsubInv = k/8; delta = 4/k; label = 'delta=4H0';
        end
        cases(end+1) = makeCase('5.8', 'impedance', 'lod', k, 'k-derived', NaN, k, HsubInv, delta, 2, ...
            label, paper58Lod(i, j)); %#ok<AGROW>
        cases(end+1) = makeCase('5.8', 'impedance', 'p1', k, 'k-derived', NaN, k, HsubInv, delta, 2, ...
            label, paper58P1{i, j}); %#ok<AGROW>
    end
end

paper59 = [9, 11, 12, 13; 8, 12, 14, 15; 8, 19, 22, 26];
ksSmallOverlap = [40, 80, 120];
for i = 1:numel(ksSmallOverlap)
    k = ksSmallOverlap(i);
    for j = 1:4
        switch j
            case 1
                delta = 1/k; label = 'delta=H0';
            case 2
                delta = '4h'; label = 'delta=4h';
            case 3
                delta = '2h'; label = 'delta=2h';
            otherwise
                delta = '1h'; label = 'delta=h';
        end
        cases(end+1) = makeCase('5.9', 'dirichlet', 'lod', k, 'fixed-kmax', NaN, k, k/2, delta, ceil(log2(k)), ...
            label, paper59(i, j)); %#ok<AGROW>
    end
end

for i = 1:numel(cases)
    cases(i).cfgCh = cfg.Ch;
end
end


function c = baseCase()
c = struct('tableId', '', 'variant', '', 'coarseType', '', 'kappa', NaN, ...
    'hMode', '', 'hInvRequested', NaN, 'targetHInv', NaN, 'targetHsubInv', NaN, ...
    'targetDelta', NaN, 'm', NaN, 'parameter', '', 'paperValue', '', ...
    'hInv', NaN, 'HInv', NaN, 'HsubInv', NaN, 'deltaSteps', NaN, ...
    'h', NaN, 'H', NaN, 'Hsub', NaN, 'delta', NaN, 'grid', [NaN, NaN], ...
    'cfgCh', NaN);
end


function c = makeCase(tableId, variant, coarseType, k, hMode, hInv, HInv, HsubInv, delta, m, parameter, paperValue)
c = baseCase();
c.tableId = tableId;
c.variant = variant;
c.coarseType = coarseType;
c.kappa = k;
c.hMode = hMode;
c.hInvRequested = hInv;
c.targetHInv = HInv;
c.targetHsubInv = HsubInv;
c.targetDelta = delta;
c.m = m;
c.parameter = parameter;
if isnumeric(paperValue)
    c.paperValue = sprintf('%g', paperValue);
else
    c.paperValue = paperValue;
end
end


function c = finalizeCase(c, cfg)
switch c.hMode
    case 'fixed'
        nh = c.hInvRequested;
    case 'fixed-kmax'
        raw = ceil(cfg.Ch * 120^(3/2));
        nh = alignFineInv(raw, [40, 60, 80, 120]);
    otherwise
        raw = ceil(cfg.Ch * c.kappa^(3/2));
        nh = alignFineInv(raw, [c.targetHInv, c.targetHsubInv]);
end

c.hInv = nh;
c.HInv = compatibleDivisor(nh, c.targetHInv);
c.HsubInv = compatibleDivisor(nh, c.targetHsubInv);

if ischar(c.targetDelta) || isstring(c.targetDelta)
    token = char(c.targetDelta);
    switch token
        case '4h'
            c.deltaSteps = 4;
        case '2h'
            c.deltaSteps = 2;
        case '1h'
            c.deltaSteps = 1;
        otherwise
            error('verify_lxzz25:delta', 'Unknown delta token "%s".', token);
    end
else
    c.deltaSteps = max(1, round(c.targetDelta * nh));
end

c.h = 1 / c.hInv;
c.H = 1 / c.HInv;
c.Hsub = 1 / c.HsubInv;
c.delta = c.deltaSteps / c.hInv;
c.grid = [round(2 * c.HsubInv), round(2 * c.HsubInv)];
end


function est = estimateCase(c, cfg)
N = (c.hInv + 1)^2;
NT = 2 * c.hInv^2;
nSub = (c.grid(1) + 1) * (c.grid(2) + 1);
localWidth = min(1, c.Hsub);
nl = max(4, (ceil(localWidth * c.hInv) + 1)^2);
globalBytes = 112 * N;
topologyBytes = 48 * NT;
gmresBytes = 16 * N * (cfg.gmresBasisLength + 3);
if strcmpi(cfg.solverMode, 'direct')
    localBytes = nSub * 336 * nl + 16 * cfg.luFillConstant * nl * log2(max(nl, 2));
else
    localBytes = nSub * (336 * nl + 16 * cfg.luFillConstant * nl * log2(max(nl, 2)));
end

if strcmpi(c.coarseType, 'lod')
    patchWidth = min(1, (2 * c.m + 1) * c.H);
    nPatch = max(4, (ceil(patchWidth * c.hInv) + 1)^2);
    nCoarseElem = 2 * c.HInv^2;
    lodBasisBytes = 2 * nCoarseElem * min(N, nPatch) * 3 * 16;
    lodPeakBytes = 2 * 16 * cfg.luFillConstant * nPatch * log2(max(nPatch, 2));
else
    nPatch = 0;
    nCoarseElem = 2 * c.HInv^2;
    lodBasisBytes = 0;
    lodPeakBytes = 0;
end
coarseBytes = 112 * (c.HInv + 1)^2;
totalBytes = globalBytes + topologyBytes + gmresBytes + localBytes + lodBasisBytes + lodPeakBytes + coarseBytes;

est = struct('N', N, 'NT', NT, 'nSub', nSub, 'nl', nl, ...
    'nPatch', nPatch, 'nCoarse', (c.HInv + 1)^2, 'nCoarseElem', nCoarseElem, ...
    'totalGB', totalBytes / 2^30, 'globalGB', globalBytes / 2^30, ...
    'localGB', localBytes / 2^30, 'lodGB', (lodBasisBytes + lodPeakBytes) / 2^30);
end


function out = runOneCase(c, cfg)
timerAll = tic;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], c.H);
parts = coarseHatPartition2D(node, elem, bdFlag, c.Hsub / 2);

lodOpts = struct('oversampling', c.m, 'solveCoarse', false, ...
    'solverMode', cfg.lodSolverMode, 'useParfor', cfg.useParfor);
preOpts = struct('variant', c.variant, 'coarseType', c.coarseType, ...
    'lodOptions', lodOpts, 'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'adjointType', cfg.adjointType);
pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, c.kappa, parts, nodeH, elemH, bdH, preOpts);

b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.kappa, 1);
applyTime = 0;
    function y = applyTimed(r)
        if toc(timerAll) > cfg.timeLimitS
            error('verify_lxzz25:timeLimit', ...
                'Experiment exceeded %.1f seconds.', cfg.timeLimitS);
        end
        t = tic;
        y = pre.applyResidual(r);
        applyTime = applyTime + toc(t);
    end

[~, flag, relres, iter, resvec] = gmres(pre.A, b, cfg.restart, cfg.tol, cfg.maxit, @applyTimed);

out = struct();
out.flag = flag;
out.finalRelres = relres;
out.gmresIter = gmresIterationCount(iter, cfg.restart);
out.resvecLength = numel(resvec);
out.localDofMin = pre.local.localDofMin;
out.localDofMax = pre.local.localDofMax;
out.localDofMean = pre.local.localDofMean;
out.lodSetupS = pre.timing.coarseSetup;
out.localSetupS = pre.timing.localSetup;
out.applyS = applyTime;
out.totalS = toc(timerAll);
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


function n = gmresIterationCount(iter, restart)
if isempty(iter)
    n = NaN;
elseif numel(iter) == 1
    n = iter;
elseif isempty(restart)
    n = iter(2);
else
    n = (iter(1) - 1) * restart + iter(2);
end
end


function r = emptyResult()
r = struct('tableId', '', 'variant', '', 'coarseType', '', 'kappa', NaN, ...
    'parameter', '', 'paperValue', '', 'hInv', NaN, 'HInv', NaN, 'HsubInv', NaN, ...
    'h', NaN, 'H', NaN, 'Hsub', NaN, 'delta', NaN, 'm', NaN, ...
    'ndof', NaN, 'nt', NaN, 'ncoarse', NaN, 'ncoarseElem', NaN, 'nsub', NaN, ...
    'localDofEst', NaN, 'patchDofEst', NaN, 'estimateGB', NaN, ...
    'localDofMin', NaN, 'localDofMax', NaN, 'localDofMean', NaN, ...
    'gmresIter', NaN, 'flag', NaN, 'finalRelres', NaN, 'resvecLength', NaN, ...
    'lodSetupS', NaN, 'localSetupS', NaN, 'applyS', NaN, 'totalS', NaN, ...
    'status', '', 'notes', '');
end


function r = copyCaseFields(r, c, est)
r.tableId = c.tableId;
r.variant = c.variant;
r.coarseType = c.coarseType;
r.kappa = c.kappa;
r.parameter = c.parameter;
r.paperValue = c.paperValue;
r.hInv = c.hInv;
r.HInv = c.HInv;
r.HsubInv = c.HsubInv;
r.h = c.h;
r.H = c.H;
r.Hsub = c.Hsub;
r.delta = c.delta;
r.m = c.m;
r.ndof = est.N;
r.nt = est.NT;
r.ncoarse = est.nCoarse;
r.ncoarseElem = est.nCoarseElem;
r.nsub = est.nSub;
r.localDofEst = est.nl;
r.patchDofEst = est.nPatch;
r.estimateGB = est.totalGB;
end


function r = mergeRunOutput(r, out)
names = fieldnames(out);
for i = 1:numel(names)
    r.(names{i}) = out.(names{i});
end
end


function writeOutputs(outDir, results, cfg)
csvPath = fullfile(outDir, 'lxzz25_article_results.csv');
mdPath = fullfile(outDir, 'lxzz25_article_results.md');
writeCsv(csvPath, results);
writeMarkdown(mdPath, results, cfg);
end


function writeCsv(path, results)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
names = fieldnames(results);
csvNames = csvHeaderNames(names);
fprintf(fid, '%s\n', strjoin(csvNames.', ','));
for i = 1:numel(results)
    vals = cell(1, numel(names));
    for j = 1:numel(names)
        v = results(i).(names{j});
        if isnumeric(v)
            vals{j} = sprintf('%.16g', v);
        else
            vals{j} = csvEscape(v);
        end
    end
    fprintf(fid, '%s\n', strjoin(vals, ','));
end
end


function csvNames = csvHeaderNames(names)
csvNames = names;
for i = 1:numel(csvNames)
    switch csvNames{i}
        case 'hInv'
            csvNames{i} = 'fine_inv';
        case 'HInv'
            csvNames{i} = 'coarse_inv';
        case 'HsubInv'
            csvNames{i} = 'subdomain_inv';
        case 'h'
            csvNames{i} = 'fine_h';
        case 'H'
            csvNames{i} = 'coarse_h';
        case 'Hsub'
            csvNames{i} = 'subdomain_h';
    end
end
end


function s = csvEscape(v)
s = char(v);
if contains(s, ',') || contains(s, '"') || contains(s, newline)
    s = ['"', strrep(s, '"', '""'), '"'];
end
end


function writeMarkdown(path, results, cfg)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Lu--Xu--Zheng--Zou (2025), Section 5 Tables 5.1--5.9.\n');
fprintf(fid, 'Created: 2026-06-01\n');
fprintf(fid, 'Updated: 2026-06-01\n');
fprintf(fid, 'Verification entry point: `verify/verify_lxzz25_article_experiments.m`\n');
fprintf(fid, 'Main utilities: `twoLevelHybridSchwarzHelmholtz2D`, `buildLODHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# LXZZ25 Article Experiment Results\n\n');
fprintf(fid, 'Memory limit: %.1f GB. Per-experiment time limit: %.0f seconds. Adjoint type: `%s`. Interactive run cap: `N <= %d` and coarse elements `<= %d`; set `LXZZ25_RUN_ALL_PERMITTED=1` to run all memory-permitted rows.\n\n', ...
    cfg.memoryLimitGB, cfg.timeLimitS, cfg.adjointType, cfg.maxRunDof, cfg.maxRunCoarseElem);
fprintf(fid, '| table | variant | coarse | parameter | kappa | h | H | Hsub | delta | m | N | estimate GB | paper | repo it | relres | status | notes |\n');
fprintf(fid, '|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    repoIt = '-';
    relres = '-';
    if ~isnan(r.gmresIter), repoIt = sprintf('%g', r.gmresIter); end
    if ~isnan(r.finalRelres), relres = sprintf('%.3e', r.finalRelres); end
    fprintf(fid, '| %s | %s | %s | %s | %.0f | 1/%d | 1/%d | 1/%d | %.4g | %d | %d | %.2f | %s | %s | %s | %s | %s |\n', ...
        r.tableId, r.variant, r.coarseType, r.parameter, r.kappa, r.hInv, r.HInv, r.HsubInv, ...
        r.delta, r.m, r.ndof, r.estimateGB, r.paperValue, repoIt, relres, r.status, markdownEscape(r.notes));
end
end


function s = markdownEscape(v)
s = char(v);
s = strrep(s, '|', '\\|');
s = strrep(s, newline, ' ');
end


function ensureParallelPool()
try
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', feature('numcores'));
    end
catch ME
    warning('verify_lxzz25:parpool', 'Could not start parpool: %s', ME.message);
end
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
