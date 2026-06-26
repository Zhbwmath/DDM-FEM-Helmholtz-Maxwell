% VERIFY_HL25_HULI_BETA_ADAPTIVE  Sweep Hu-Li subdomain beta with LXZZ locals.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'HL25_Helmholtz_harmonic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg = localConfig();
nCase = numel(cfg.kValues) * numel(cfg.betaValues);
rows = repmat(emptyRow(), nCase, 1);
rowCount = 0;
fprintf('========== Hu-Li beta adaptive sweep ==========\n');
fprintf('k values: %s, beta values: %s\n', mat2str(cfg.kValues), mat2str(cfg.betaValues));

for ik = 1:numel(cfg.kValues)
    for ib = 1:numel(cfg.betaValues)
        c = finalizeCase(cfg.kValues(ik), cfg.betaValues(ib), cfg);
        est = estimateCase(c, cfg);
        row = baseRow(c, est, cfg);
        fprintf('k=%g beta=%.3g: N=%d coarseHInv=%d raw=%d est %.2f GB ... ', ...
            c.k, c.beta, est.N, c.coarseHInv, est.rawCoarseEstimate, est.totalGB);
        if est.N > cfg.maxRunDof
            row.status = 'estimated_dof_cap';
            row.notes = sprintf('Above adaptive run cap N=%d.', cfg.maxRunDof);
            fprintf('ESTIMATED DOF\n');
        elseif est.rawCoarseEstimate > cfg.maxCoarseEstimate
            row.status = 'estimated_coarse_cap';
            row.notes = sprintf('Above adaptive coarse cap %d.', cfg.maxCoarseEstimate);
            fprintf('ESTIMATED COARSE\n');
        elseif ~cfg.runEnabled
            row.status = 'estimated_only';
            row.notes = 'Execution disabled by HL25_HULI_BETA_RUN=0.';
            fprintf('ESTIMATED\n');
        else
            try
                row = runCase(row, c, cfg);
                fprintf('RAN: %d iterations, relres %.3e\n', ...
                    row.gmresIterations, row.finalRelres);
            catch ME
                row.status = 'failed';
                row.notes = sprintf('%s: %s', ME.identifier, ME.message);
                fprintf('FAILED: %s\n', ME.message);
            end
        end
        rowCount = rowCount + 1;
        rows(rowCount) = row;
        writeOutputs(outDir, rows(1:rowCount), cfg);
    end
end
writeOutputs(outDir, rows(1:rowCount), cfg);
fprintf('CSV: %s\n', fullfile(outDir, [cfg.outputStem, '.csv']));
fprintf('Markdown: %s\n', fullfile(outDir, [cfg.outputStem, '.md']));
fprintf('========== Hu-Li beta adaptive sweep complete ==========\n');


function cfg = localConfig()
cfg.kValues = envVector('HL25_HULI_BETA_KVALUES', [32, 64]);
cfg.betaValues = envVector('HL25_HULI_BETA_VALUES', [0.4, 0.5, 0.6, 0.7]);
cfg.degree = envNumber('HL25_HULI_BETA_DEGREE', 1);
cfg.Ch = envNumber('HL25_HULI_BETA_CH', 1);
cfg.variant = envString('HL25_HULI_BETA_VARIANT', 'dirichlet');
cfg.coarseType = envString('HL25_HULI_BETA_COARSE', 'economic');
cfg.tol = envNumber('HL25_HULI_BETA_TOL', 1e-6);
cfg.maxit = envNumber('HL25_HULI_BETA_MAXIT', 100);
cfg.runEnabled = logical(envNumber('HL25_HULI_BETA_RUN', 1));
cfg.maxRunDof = envNumber('HL25_HULI_BETA_MAX_RUN_DOF', 80000);
cfg.maxCoarseEstimate = envNumber('HL25_HULI_BETA_MAX_COARSE', 30000);
cfg.solverMode = envString('HL25_HULI_BETA_SOLVER_MODE', 'adaptive');
cfg.useParfor = logical(envNumber('HL25_HULI_BETA_PARFOR', 0));
cfg.localLuFillConstant = envNumber('HL25_HULI_BETA_LU_FILL', 20);
cfg.outputStem = envString('HL25_HULI_BETA_OUTPUT_STEM', ...
    'huli_beta_adaptive_results');
end


function c = finalizeCase(k, beta, cfg)
c = struct();
c.k = k;
c.beta = beta;
c.degree = cfg.degree;
c.nu = max(1, round(k^(1 - beta)));
c.rho = 0.5 * k^((beta - 1) / 2);
c.localHInv = max(1, round(k));
c.coarseHInv = max(1, round(k^beta));
rawFine = ceil(cfg.Ch * k^(3/2));
c.hInv = alignFineInv(rawFine, lcm(c.localHInv, c.coarseHInv));
c.h = 1 / c.hInv;
c.localOverlap = max(1, round((1 / c.localHInv) * c.hInv)) / c.hInv;
c.coarseOverlap = max(1, round((0.25 / c.coarseHInv) * c.hInv)) / c.hInv;
c.nCoarseSubdomains = c.coarseHInv^2;
c.nLocalSubdomains = c.localHInv^2;
end


function est = estimateCase(c, cfg)
if c.degree == 1
    N = (c.hInv + 1)^2;
else
    N = (2 * c.hInv + 1)^2;
end
NT = 2 * c.hInv^2;
supportSteps = max(1, ceil((1 / c.coarseHInv + 2 * c.coarseOverlap) * c.hInv));
if c.degree == 1
    coarseLocalDof = (supportSteps + 1)^2;
else
    coarseLocalDof = (2 * supportSteps + 1)^2;
end
if strcmpi(cfg.coarseType, 'economic')
    rawCoarse = c.nCoarseSubdomains * (2 * c.nu);
else
    rawCoarse = c.nCoarseSubdomains * max(4, 4 * supportSteps);
end
globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.maxit + 3);
localBytes = c.nCoarseSubdomains * 336 * coarseLocalDof;
localLuBytes = c.nCoarseSubdomains * 16 * cfg.localLuFillConstant * ...
    coarseLocalDof * log2(max(coarseLocalDof, 2));
coarseBytes = rawCoarse * 80 * 16;
coarseLuBytes = 16 * cfg.localLuFillConstant * rawCoarse * ...
    log2(max(rawCoarse, 2));
est = struct('N', N, 'NT', NT, 'coarseSupportSteps', supportSteps, ...
    'coarseLocalDofEstimate', coarseLocalDof, ...
    'rawCoarseEstimate', rawCoarse, ...
    'totalGB', (globalBytes + gmresBytes + localBytes + localLuBytes + ...
    coarseBytes + coarseLuBytes) / 2^30);
end


function row = runCase(row, c, cfg)
tAll = tic;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
bbox = [0, 1, 0, 1];
coarseParts = rectangularParts(node, elem, bbox, [c.coarseHInv, c.coarseHInv], ...
    c.coarseOverlap);
coarseParts = linearPartitionOfUnity2D(coarseParts, bbox, ...
    [c.coarseHInv, c.coarseHInv], c.coarseOverlap);
localParts = rectangularParts(node, elem, bbox, [c.localHInv, c.localHInv], ...
    c.localOverlap);
localParts = linearPartitionOfUnity2D(localParts, bbox, ...
    [c.localHInv, c.localHInv], c.localOverlap);
pde = helmholtzPDE(c.k, 'epsilon', 0);
builderOpts = struct('degree', c.degree, 'coarseType', cfg.coarseType, ...
    'rho', c.rho, 'nu', c.nu, 'kappaRef', c.k, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'rankMethod', 'none', 'coarseSolverMode', 'lu', ...
    'localLuFillConstant', cfg.localLuFillConstant, ...
    'cacheEnergySolver', true, 'cacheEnergyAdjoint', false);
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, coarseParts, builderOpts);
preOpts = struct('fineSpace', method.fineSpace, ...
    'coarseSpace', method.coarseSpace, 'variant', cfg.variant, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'adjointType', 'energy', ...
    'localLuFillConstant', cfg.localLuFillConstant);
pre = twoLevelHybridSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, localParts, [], [], [], preOpts);
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.k, c.degree);
tSolve = tic;
[~, flag, relres, iter, resvec] = gmres( ...
    pre.A, b, [], cfg.tol, cfg.maxit, @pre.applyResidual);
row.flag = flag;
row.finalRelres = relres;
row.gmresIterations = gmresIterationCount(iter);
row.resvecLength = numel(resvec);
row.coarseDimension = method.stats.coarseDimension;
row.rawCoarseDimension = method.stats.rawCoarseDimension;
row.setupS = method.timing.totalSetup + pre.timing.localSetup;
row.solveS = toc(tSolve);
row.totalS = toc(tAll);
row.status = ternary(flag == 0, 'ran', 'ran_not_converged');
row.notes = '';
end


function row = emptyRow()
row = struct('k', NaN, 'beta', NaN, 'degree', NaN, 'coarseType', '', ...
    'variant', '', 'nu', NaN, 'rho', NaN, 'fineHInv', NaN, ...
    'localHInv', NaN, 'coarseHInv', NaN, 'localOverlap', NaN, ...
    'coarseOverlap', NaN, 'ndof', NaN, 'nLocalSubdomains', NaN, ...
    'nCoarseSubdomains', NaN, 'coarseLocalDofEstimate', NaN, ...
    'rawCoarseEstimate', NaN, 'estimateGB', NaN, ...
    'coarseDimension', NaN, 'rawCoarseDimension', NaN, ...
    'gmresIterations', NaN, 'flag', NaN, 'finalRelres', NaN, ...
    'resvecLength', NaN, 'setupS', NaN, 'solveS', NaN, ...
    'totalS', NaN, 'status', '', 'notes', '');
end


function row = baseRow(c, est, cfg)
row = emptyRow();
row.k = c.k;
row.beta = c.beta;
row.degree = c.degree;
row.coarseType = cfg.coarseType;
row.variant = cfg.variant;
row.nu = c.nu;
row.rho = c.rho;
row.fineHInv = c.hInv;
row.localHInv = c.localHInv;
row.coarseHInv = c.coarseHInv;
row.localOverlap = c.localOverlap;
row.coarseOverlap = c.coarseOverlap;
row.ndof = est.N;
row.nLocalSubdomains = c.nLocalSubdomains;
row.nCoarseSubdomains = c.nCoarseSubdomains;
row.coarseLocalDofEstimate = est.coarseLocalDofEstimate;
row.rawCoarseEstimate = est.rawCoarseEstimate;
row.estimateGB = est.totalGB;
row.status = 'pending';
end


function writeOutputs(outDir, rows, cfg)
csvPath = fullfile(outDir, [cfg.outputStem, '.csv']);
mdPath = fullfile(outDir, [cfg.outputStem, '.md']);
names = fieldnames(emptyRow());
fid = fopen(csvPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(names, ','));
for i = 1:numel(rows)
    vals = cell(1, numel(names));
    for j = 1:numel(names)
        v = rows(i).(names{j});
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
fprintf(fid, 'Reproduction target: Hu--Li subdomain-size adaptive sweep with LXZZ local solver settings.\n');
fprintf(fid, 'Created: 2026-06-26\nUpdated: 2026-06-26\n');
fprintf(fid, 'Verification entry point: `verify/verify_hl25_huli_beta_adaptive.m`\n');
fprintf(fid, 'Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, structured rectangular element binning, MATLAB `gmres`\n\n');
fprintf(fid, '# Hu--Li Beta Adaptive Sweep\n\n');
fprintf(fid, 'Hu--Li coarse spaces use $d\\approx\\kappa^{-\\beta}$ and overlap $d/4$; LXZZ local solvers keep $H=1/\\kappa$ and overlap $H$. Rows above the adaptive caps are estimates only.\n\n');
fprintf(fid, '| k | beta | h^{-1} | coarse H^{-1} | coarse sub | raw coarse | coarse dim | est GB | iter | relres | status | notes |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(rows)
    r = rows(i);
    fprintf(fid, '| %g | %.3g | %g | %g | %g | %g | %s | %.2f | %s | %s | %s | %s |\n', ...
        r.k, r.beta, r.fineHInv, r.coarseHInv, r.nCoarseSubdomains, ...
        r.rawCoarseEstimate, mdNumber(firstFinite(r.coarseDimension, ...
        r.rawCoarseEstimate)), r.estimateGB, mdNumber(r.gmresIterations), ...
        mdNumber(r.finalRelres), r.status, mdEscape(r.notes));
end
clear cleanup
end


function parts = rectangularParts(node, elem, bbox, gridSize, overlap)
nx = gridSize(1); ny = gridSize(2);
nTotal = nx * ny;
xmin = bbox(1); xmax = bbox(2); ymin = bbox(3); ymax = bbox(4);
Hx = (xmax - xmin) / nx; Hy = (ymax - ymin) / ny;
vertexElem = elem(:, 1:3);
xC = mean(reshape(node(vertexElem, 1), size(vertexElem)), 2);
yC = mean(reshape(node(vertexElem, 2), size(vertexElem)), 2);
ix = min(nx, max(1, floor((xC - xmin) / Hx) + 1));
iy = min(ny, max(1, floor((yC - ymin) / Hy) + 1));
binId = (iy - 1) * nx + ix;
elemByBin = accumarray(binId, (1:size(vertexElem, 1))', [nTotal, 1], @(v) {v});
elemByBin(cellfun(@isempty, elemByBin)) = {zeros(0, 1)};
parts = repmat(struct('elemIdx', []), nTotal, 1);
tol = max([xmax - xmin, ymax - ymin, 1]) * 1e-12;
for j = 1:ny
    for i = 1:nx
        s = (j - 1) * nx + i;
        xL = xmin + (i - 1) * Hx; xR = xmin + i * Hx;
        yB = ymin + (j - 1) * Hy; yT = ymin + j * Hy;
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
        keep = xC(candidates) >= xExtL - tol & xC(candidates) <= xExtR + tol & ...
            yC(candidates) >= yExtB - tol & yC(candidates) <= yExtT + tol;
        parts(s).elemIdx = candidates(keep);
    end
end
end


function hInv = alignFineInv(raw, divisor)
hInv = max(1, ceil(raw));
if divisor > 0
    hInv = ceil(hInv / divisor) * divisor;
end
end


function n = gmresIterationCount(iter)
if numel(iter) == 2
    n = iter(2);
else
    n = iter;
end
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


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


function v = firstFinite(a, b)
if isfinite(a), v = a; else, v = b; end
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
