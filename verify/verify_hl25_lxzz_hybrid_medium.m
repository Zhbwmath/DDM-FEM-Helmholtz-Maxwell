% VERIFY_HL25_LXZZ_HYBRID_MEDIUM  Hu-Li coarse spaces in LXZZ settings.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'HL25_Helmholtz_harmonic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg = localConfig();
cases = localCases(cfg);
if cfg.resume
    results = loadResults(outDir, cfg);
else
    results = repmat(emptyRow(), 0, 1);
end

fprintf('========== Hu-Li coarse space with LXZZ hybrid two-level DDM ==========\n\n');
fprintf('k values: %s\n', mat2str(cfg.kValues));
fprintf(['degree P%d, C_h=%g, beta=%g, rectangular overlap factor=%g, ' ...
    'workers=%g, max run dof=%g, max coarse estimate=%g\n\n'], ...
    cfg.degree, cfg.Ch, cfg.beta, cfg.overlapFactor, ...
    effectiveWorkerCount(cfg), cfg.maxRunDof, cfg.maxCoarseEstimate);

for i = 1:numel(cases)
    c = finalizeCase(cases(i), cfg);
    if cfg.resume && caseComplete(results, c)
        fprintf('[%02d/%02d] k=%g coarse=%s variant=%s: checkpoint complete, skipping.\n', ...
            i, numel(cases), c.k, c.coarseType, c.variant);
        continue;
    end
    est = estimateCase(c, cfg);
    fprintf('[%02d/%02d] k=%g coarse=%s variant=%s: N=%d, coarse est=%d, %.2f GB ... ', ...
        i, numel(cases), c.k, c.coarseType, c.variant, est.N, ...
        est.rawCoarseEstimate, est.totalGB);

    row = baseRow(c, est);
    if est.totalGB > cfg.hardMemoryGB
        row.status = 'blocked_memory_gt_hard_limit';
        row.notes = sprintf('Estimate %.2f GB exceeds hard limit %.2f GB.', ...
            est.totalGB, cfg.hardMemoryGB);
        fprintf('BLOCKED\n');
    elseif est.totalGB > cfg.permissionMemoryGB && ~cfg.allowPermissionRows
        row.status = 'requires_permission_gt_200gb';
        row.notes = sprintf('Estimate %.2f GB requires permission gate %.2f GB.', ...
            est.totalGB, cfg.permissionMemoryGB);
        fprintf('PERMISSION REQUIRED\n');
    elseif est.N > cfg.maxRunDof
        row.status = 'queued_runtime_cap';
        row.notes = sprintf('Above medium run cap N=%d.', cfg.maxRunDof);
        fprintf('QUEUED DOF\n');
    elseif est.rawCoarseEstimate > cfg.maxCoarseEstimate
        row.status = 'queued_coarse_cap';
        row.notes = sprintf('Above medium coarse cap %d.', cfg.maxCoarseEstimate);
        fprintf('QUEUED COARSE\n');
    elseif ~cfg.runEnabled
        row.status = 'estimated_only';
        row.notes = 'Execution disabled by HL25_LXZZ_MEDIUM_RUN=0.';
        fprintf('ESTIMATED\n');
    else
        try
            row = runCase(row, c, est, cfg);
            fprintf('RAN: %d iterations, relres %.3e\n', ...
                row.gmresIterations, row.finalRelres);
        catch ME
            row.status = 'failed';
            row.notes = sprintf('%s: %s', ME.identifier, ME.message);
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    results = replaceRow(results, row);
    writeOutputs(outDir, results, cfg);
end

writeOutputs(outDir, results, cfg);
fprintf('\nCSV: %s\n', fullfile(outDir, [cfg.outputStem, '.csv']));
fprintf('Markdown: %s\n', fullfile(outDir, [cfg.outputStem, '.md']));
fprintf('========== Hu-Li/LXZZ medium verification complete ==========\n');


function cfg = localConfig()
cfg.kValues = envVector('HL25_LXZZ_MEDIUM_KVALUES', [16, 32, 64, 128]);
cfg.coarseTypes = envList('HL25_LXZZ_MEDIUM_COARSE', {'economic', 'spectral'});
cfg.variants = envList('HL25_LXZZ_MEDIUM_VARIANTS', {'dirichlet', 'impedance'});
cfg.degree = envNumber('HL25_LXZZ_MEDIUM_DEGREE', 1);
cfg.beta = envNumber('HL25_LXZZ_MEDIUM_BETA', 0.6);
cfg.Ch = envNumber('HL25_LXZZ_MEDIUM_CH', 1);
cfg.overlapFactor = envNumber('HL25_LXZZ_MEDIUM_OVERLAP_FACTOR', 1);
cfg.tol = envNumber('HL25_LXZZ_MEDIUM_TOL', 1e-6);
cfg.maxit = envNumber('HL25_LXZZ_MEDIUM_MAXIT', 100);
cfg.restart = [];
cfg.runEnabled = logical(envNumber('HL25_LXZZ_MEDIUM_RUN', 1));
cfg.resume = logical(envNumber('HL25_LXZZ_MEDIUM_RESUME', 1));
cfg.maxRunDof = envNumber('HL25_LXZZ_MEDIUM_MAX_RUN_DOF', 50000);
cfg.maxCoarseEstimate = envNumber('HL25_LXZZ_MEDIUM_MAX_COARSE', 20000);
cfg.permissionMemoryGB = envNumber('HL25_LXZZ_MEDIUM_PERMISSION_GB', 200);
cfg.hardMemoryGB = envNumber('HL25_LXZZ_MEDIUM_HARD_GB', 500);
cfg.allowPermissionRows = logical(envNumber('HL25_LXZZ_MEDIUM_ALLOW_GT_200', 0));
cfg.solverMode = envString('HL25_LXZZ_MEDIUM_SOLVER_MODE', 'adaptive');
cfg.coarseSolverMode = envString('HL25_LXZZ_MEDIUM_COARSE_SOLVER', 'lu');
cfg.adjointType = envString('HL25_LXZZ_MEDIUM_ADJOINT', 'energy');
cfg.rankMethod = envString('HL25_LXZZ_MEDIUM_RANK_METHOD', 'none');
cfg.useParfor = logical(envNumber('HL25_LXZZ_MEDIUM_PARFOR', 0));
cfg.workerCount = envNumber('HL25_LXZZ_MEDIUM_WORKERS', 0);
cfg.localLuFillConstant = envNumber('HL25_LXZZ_MEDIUM_LU_FILL', 20);
cfg.coarseLuFillConstant = envNumber('HL25_LXZZ_MEDIUM_COARSE_LU_FILL', 20);
cfg.gmresBasisLength = envNumber('HL25_LXZZ_MEDIUM_GMRES_BASIS', 103);
cfg.outputStem = envString('HL25_LXZZ_MEDIUM_OUTPUT_STEM', ...
    'lxzz_hl_coarse_medium_results');
end


function cases = localCases(cfg)
cases = repmat(emptyCase(), 0, 1);
for ik = 1:numel(cfg.kValues)
    for ic = 1:numel(cfg.coarseTypes)
        for iv = 1:numel(cfg.variants)
            c = emptyCase();
            c.k = cfg.kValues(ik);
            c.coarseType = cfg.coarseTypes{ic};
            c.variant = cfg.variants{iv};
            cases(end+1) = c; %#ok<AGROW>
        end
    end
end
end


function c = emptyCase()
c = struct('k', NaN, 'degree', NaN, 'coarseType', '', 'variant', '', ...
    'beta', NaN, 'rho', NaN, 'nu', NaN, 'hInv', NaN, 'h', NaN, ...
    'HInv', NaN, 'H', NaN, 'overlap', NaN, 'nSubdomains', NaN);
end


function c = finalizeCase(c, cfg)
c.degree = cfg.degree;
c.beta = cfg.beta;
c.rho = 0.5 * c.k^((c.beta - 1) / 2);
c.nu = max(1, round(c.k^(1 - c.beta)));
c.HInv = max(1, round(c.k));
rawFine = ceil(cfg.Ch * c.k^(3/2));
c.hInv = alignFineInv(rawFine, c.HInv);
c.h = 1 / c.hInv;
c.H = 1 / c.HInv;
overlapCells = max(1, round(cfg.overlapFactor * c.H * c.hInv));
c.overlap = overlapCells / c.hInv;
c.nSubdomains = c.HInv^2;
end


function est = estimateCase(c, cfg)
if c.degree == 1
    N = (c.hInv + 1)^2;
else
    N = (2 * c.hInv + 1)^2;
end
NT = 2 * c.hInv^2;
supportSteps = max(1, ceil((c.H + 2 * c.overlap) * c.hInv));
if c.degree == 1
    localDof = (supportSteps + 1)^2;
    spectralBoundaryEstimate = max(4, 4 * supportSteps);
else
    localDof = (2 * supportSteps + 1)^2;
    spectralBoundaryEstimate = max(8, 8 * supportSteps);
end
if strcmpi(c.coarseType, 'economic')
    rawCoarse = c.nSubdomains * (2 * c.nu);
    modesPerSubdomain = 2 * c.nu;
else
    rawCoarse = c.nSubdomains * spectralBoundaryEstimate;
    modesPerSubdomain = spectralBoundaryEstimate;
end
globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.gmresBasisLength + 3);
localBytes = c.nSubdomains * 336 * localDof;
localLuBytes = c.nSubdomains * 16 * cfg.localLuFillConstant * ...
    localDof * log2(max(localDof, 2));
oneLocalLuBytes = 336 * localDof + ...
    16 * cfg.localLuFillConstant * localDof * log2(max(localDof, 2));
basisBytesPerSubdomain = 16 * localDof * modesPerSubdomain;
workers = effectiveWorkerCount(cfg);
subdomainsPerWorker = ceil(c.nSubdomains / workers);
workerBytes = globalBytes + oneLocalLuBytes + ...
    subdomainsPerWorker * basisBytesPerSubdomain;
coarseBytes = rawCoarse * 80 * 16;
coarseLuBytes = 16 * cfg.coarseLuFillConstant * rawCoarse * ...
    log2(max(rawCoarse, 2));
totalBytes = globalBytes + gmresBytes + localBytes + localLuBytes + ...
    coarseBytes + coarseLuBytes;
est = struct('N', N, 'NT', NT, 'supportSteps', supportSteps, ...
    'localDofEstimate', localDof, 'rawCoarseEstimate', rawCoarse, ...
    'workerCountEstimate', workers, ...
    'estimatedPerWorkerGB', workerBytes / 2^30, ...
    'estimatedParallelGB', (totalBytes + workers * workerBytes) / 2^30, ...
    'totalGB', totalBytes / 2^30);
end


function row = runCase(row, c, est, cfg)
tAll = tic;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
parts = rectangularLxzzPartition(node, elem, c);
pde = helmholtzPDE(c.k, 'epsilon', 0);
builderOpts = struct('degree', c.degree, 'coarseType', c.coarseType, ...
    'rho', c.rho, 'nu', c.nu, 'kappaRef', c.k, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'rankMethod', cfg.rankMethod, ...
    'coarseSolverMode', cfg.coarseSolverMode, ...
    'localLuFillConstant', cfg.localLuFillConstant, ...
    'cacheEnergySolver', strcmpi(cfg.adjointType, 'energy'), ...
    'cacheEnergyAdjoint', false);
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, parts, builderOpts);
preOpts = struct('fineSpace', method.fineSpace, ...
    'coarseSpace', method.coarseSpace, 'variant', c.variant, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'adjointType', cfg.adjointType, ...
    'localLuFillConstant', cfg.localLuFillConstant);
pre = twoLevelHybridSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, parts, [], [], [], preOpts);
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.k, c.degree);
tSolve = tic;
[~, flag, relres, iter, resvec] = gmres( ...
    pre.A, b, cfg.restart, cfg.tol, cfg.maxit, @pre.applyResidual);
row.flag = flag;
row.finalRelres = relres;
row.gmresIterations = gmresIterationCount(iter, cfg.restart);
row.resvecLength = numel(resvec);
row.coarseDimension = method.stats.coarseDimension;
row.rawCoarseDimension = method.stats.rawCoarseDimension;
row.setupS = method.timing.totalSetup + pre.timing.localSetup;
row.solveS = toc(tSolve);
row.totalS = toc(tAll);
row.localDofMin = pre.local.localDofMin;
row.localDofMax = pre.local.localDofMax;
row.localDofMean = pre.local.localDofMean;
row.localMode = pre.local.solverModeEffective;
row.status = ternary(flag == 0, 'ran', 'ran_not_converged');
row.notes = '';
if est.rawCoarseEstimate ~= method.stats.rawCoarseDimension
    row.notes = sprintf('Estimated raw coarse %d; actual %d.', ...
        est.rawCoarseEstimate, method.stats.rawCoarseDimension);
end
end


function row = emptyRow()
row = struct('k', NaN, 'degree', NaN, 'coarseType', '', 'variant', '', ...
    'beta', NaN, 'rho', NaN, 'nu', NaN, ...
    'fineHInv', NaN, 'localHInv', NaN, 'overlap', NaN, ...
    'ndof', NaN, 'nSubdomains', NaN, 'supportSteps', NaN, ...
    'localDofEstimate', NaN, 'rawCoarseEstimate', NaN, ...
    'workerCount', NaN, 'estimatedPerWorkerGB', NaN, ...
    'estimatedParallelGB', NaN, ...
    'estimateGB', NaN, 'coarseDimension', NaN, ...
    'rawCoarseDimension', NaN, 'localMode', '', 'localDofMin', NaN, ...
    'localDofMax', NaN, 'localDofMean', NaN, ...
    'gmresIterations', NaN, 'flag', NaN, 'finalRelres', NaN, ...
    'resvecLength', NaN, 'setupS', NaN, 'solveS', NaN, 'totalS', NaN, ...
    'status', '', 'notes', '');
end


function row = baseRow(c, est)
row = emptyRow();
fields = fieldnames(c);
for i = 1:numel(fields)
    if isfield(row, fields{i}), row.(fields{i}) = c.(fields{i}); end
end
row.ndof = est.N;
row.fineHInv = c.hInv;
row.localHInv = c.HInv;
row.supportSteps = est.supportSteps;
row.localDofEstimate = est.localDofEstimate;
row.rawCoarseEstimate = est.rawCoarseEstimate;
row.workerCount = est.workerCountEstimate;
row.estimatedPerWorkerGB = est.estimatedPerWorkerGB;
row.estimatedParallelGB = est.estimatedParallelGB;
row.estimateGB = est.totalGB;
row.status = 'pending';
end


function results = loadResults(outDir, cfg)
path = fullfile(outDir, [cfg.outputStem, '.csv']);
if ~exist(path, 'file')
    results = repmat(emptyRow(), 0, 1);
    return;
end
t = readtable(path, 'TextType', 'string');
results = repmat(emptyRow(), height(t), 1);
for i = 1:height(t)
    r = emptyRow();
    names = fieldnames(r);
    for j = 1:numel(names)
        if ~ismember(names{j}, t.Properties.VariableNames), continue; end
        v = t.(names{j})(i);
        if isnumeric(r.(names{j}))
            r.(names{j}) = double(v);
        elseif ismissing(v)
            r.(names{j}) = '';
        else
            r.(names{j}) = char(v);
        end
    end
    results(i) = r;
end
end


function tf = caseComplete(results, c)
tf = false;
for i = 1:numel(results)
    r = results(i);
    if r.k == c.k && strcmp(r.coarseType, c.coarseType) && ...
            strcmp(r.variant, c.variant) && ...
            ismember(r.status, {'ran', 'ran_not_converged', ...
            'queued_runtime_cap', 'queued_coarse_cap', ...
            'requires_permission_gt_200gb', 'blocked_memory_gt_hard_limit'})
        tf = true;
        return;
    end
end
end


function results = replaceRow(results, row)
for i = 1:numel(results)
    r = results(i);
    if r.k == row.k && strcmp(r.coarseType, row.coarseType) && ...
            strcmp(r.variant, row.variant)
        results(i) = row;
        return;
    end
end
results(end+1) = row;
end


function writeOutputs(outDir, results, cfg)
csvPath = fullfile(outDir, [cfg.outputStem, '.csv']);
mdPath = fullfile(outDir, [cfg.outputStem, '.md']);
if isempty(results), return; end
names = fieldnames(emptyRow());
fid = fopen(csvPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(names, ','));
for i = 1:numel(results)
    vals = cell(1, numel(names));
    for j = 1:numel(names)
        v = results(i).(names{j});
        if isnumeric(v)
            vals{j} = num2str(v, 16);
        else
            vals{j} = csvEscape(v);
        end
    end
    fprintf(fid, '%s\n', strjoin(vals, ','));
end
clear cleanup

fid = fopen(mdPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Hu--Li Helmholtz-harmonic coarse spaces in LXZZ hybrid two-level DDM settings.\n');
fprintf(fid, 'Created: 2026-06-25\n');
fprintf(fid, 'Updated: 2026-06-26\n');
fprintf(fid, 'Verification entry point: `verify/verify_hl25_lxzz_hybrid_medium.m`\n');
fprintf(fid, 'Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, structured rectangular element binning, `linearPartitionOfUnity2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# Hu--Li Coarse Space With LXZZ Hybrid Two-Level DDM\n\n');
fprintf(fid, 'Settings: P%d fine space, literal wave numbers $\\kappa\\in%s$, LXZZ local partition uses overlapping rectangular subdomains with base size $H=1/\\kappa$ and overlap $\\delta=%gH$ aligned to the fine mesh. The POU is the normalized tensor-product linear weight from `linearPartitionOfUnity2D`. The fine grid is $h^{-1}=\\operatorname{align}(\\lceil C_h\\kappa^{3/2}\\rceil,\\kappa)$ with $C_h=%g$, $\\beta=%g$, $\\epsilon=0$, adjoint `%s`, GMRES tolerance %.1e and max iterations %d. Medium gates: $N\\le %g$ and raw coarse dimension estimate $\\le %g$. Worker-count estimate uses %g workers and reports optimistic per-worker payload plus client retained storage.\n\n', ...
    cfg.degree, mat2str(cfg.kValues), cfg.overlapFactor, cfg.Ch, cfg.beta, ...
    cfg.adjointType, cfg.tol, cfg.maxit, cfg.maxRunDof, ...
    cfg.maxCoarseEstimate, effectiveWorkerCount(cfg));
fprintf(fid, '| k | coarse | variant | h^{-1} | H^{-1} | overlap | N | subdomains | coarse dim | client GB | worker GB | parallel GB | iter | relres | status | notes |\n');
fprintf(fid, '|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %g | %s | %s | %g | %g | %.4g | %g | %g | %s | %.2f | %.2f | %.2f | %s | %s | %s | %s |\n', ...
        r.k, r.coarseType, r.variant, r.fineHInv, r.localHInv, ...
        r.overlap, r.ndof, r.nSubdomains, mdNumber(firstFinite(r.coarseDimension, ...
        r.rawCoarseEstimate)), r.estimateGB, r.estimatedPerWorkerGB, ...
        r.estimatedParallelGB, mdNumber(r.gmresIterations), ...
        mdNumber(r.finalRelres), r.status, mdEscape(r.notes));
end
clear cleanup
end


function hInv = alignFineInv(raw, divisor)
hInv = max(1, ceil(raw));
if divisor > 0
    hInv = ceil(hInv / divisor) * divisor;
end
end


function workers = effectiveWorkerCount(cfg)
if cfg.useParfor
    workers = cfg.workerCount;
    if ~isfinite(workers) || workers <= 0
        workers = feature('numcores');
    end
else
    workers = 1;
end
workers = max(1, round(workers));
end


function parts = rectangularLxzzPartition(node, elem, c)
gridSize = [c.HInv, c.HInv];
bbox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
parts = structuredRectangularParts(node, elem, bbox, gridSize, c.overlap);
parts = linearPartitionOfUnity2D(parts, bbox, gridSize, c.overlap);
end


function parts = structuredRectangularParts(node, elem, bbox, gridSize, overlap)
nx = gridSize(1);
ny = gridSize(2);
nTotal = nx * ny;
xmin = bbox(1); xmax = bbox(2);
ymin = bbox(3); ymax = bbox(4);
Hx = (xmax - xmin) / nx;
Hy = (ymax - ymin) / ny;
vertexElem = elem(:, 1:3);
xC = mean(reshape(node(vertexElem, 1), size(vertexElem)), 2);
yC = mean(reshape(node(vertexElem, 2), size(vertexElem)), 2);
ix = min(nx, max(1, floor((xC - xmin) / Hx) + 1));
iy = min(ny, max(1, floor((yC - ymin) / Hy) + 1));
binId = (iy - 1) * nx + ix;
elemByBin = accumarray(binId, (1:size(vertexElem, 1))', ...
    [nTotal, 1], @(v) {v});
emptyBins = cellfun(@isempty, elemByBin);
elemByBin(emptyBins) = {zeros(0, 1)};

parts = repmat(struct('elemIdx', []), nTotal, 1);
tol = max([xmax - xmin, ymax - ymin, 1]) * 1e-12;
for j = 1:ny
    for i = 1:nx
        s = (j - 1) * nx + i;
        xL = xmin + (i - 1) * Hx;
        xR = xmin + i * Hx;
        yB = ymin + (j - 1) * Hy;
        yT = ymin + j * Hy;
        xExtL = max(xmin, xL - overlap);
        xExtR = min(xmax, xR + overlap);
        yExtB = max(ymin, yB - overlap);
        yExtT = min(ymax, yT + overlap);
        i0 = min(nx, max(1, floor((xExtL - xmin) / Hx) + 1));
        i1 = min(nx, max(1, floor((xExtR - xmin) / Hx) + 1));
        j0 = min(ny, max(1, floor((yExtB - ymin) / Hy) + 1));
        j1 = min(ny, max(1, floor((yExtT - ymin) / Hy) + 1));
        ids = zeros((i1 - i0 + 1) * (j1 - j0 + 1), 1);
        count = 0;
        for jj = j0:j1
            for ii = i0:i1
                count = count + 1;
                ids(count) = (jj - 1) * nx + ii;
            end
        end
        candidates = vertcat(elemByBin{ids});
        keep = xC(candidates) >= xExtL - tol & ...
            xC(candidates) <= xExtR + tol & ...
            yC(candidates) >= yExtB - tol & ...
            yC(candidates) <= yExtT + tol;
        parts(s).elemIdx = candidates(keep);
    end
end
end


function n = gmresIterationCount(iter, restart)
if isempty(restart)
    if numel(iter) == 2
        n = iter(2);
    else
        n = iter;
    end
else
    n = (iter(1) - 1) * restart + iter(2);
end
end


function v = firstFinite(a, b)
if isfinite(a), v = a; else, v = b; end
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


function s = envString(name, defaultValue)
s = getenv(name);
if isempty(s), s = defaultValue; end
end


function n = envNumber(name, defaultValue)
txt = getenv(name);
if isempty(txt)
    n = defaultValue;
else
    n = str2double(txt);
    if isnan(n), n = defaultValue; end
end
end


function v = envVector(name, defaultValue)
txt = getenv(name);
if isempty(txt)
    v = defaultValue;
    return;
end
parts = split(string(txt), ',');
v = zeros(1, numel(parts));
for i = 1:numel(parts)
    v(i) = str2double(parts(i));
end
v = v(isfinite(v));
if isempty(v), v = defaultValue; end
end


function v = envList(name, defaultValue)
txt = getenv(name);
if isempty(txt)
    v = defaultValue;
    return;
end
parts = split(string(txt), ',');
v = cell(1, numel(parts));
for i = 1:numel(parts)
    v{i} = char(strtrim(parts(i)));
end
v = v(~cellfun(@isempty, v));
if isempty(v), v = defaultValue; end
end


function s = csvEscape(v)
s = char(string(v));
if contains(s, ',') || contains(s, '"')
    s = ['"', strrep(s, '"', '""'), '"'];
end
end


function s = mdNumber(v)
if isnumeric(v) && isfinite(v)
    s = sprintf('%.4g', v);
else
    s = '-';
end
end


function s = mdEscape(v)
s = char(string(v));
s = strrep(s, '|', '/');
end
