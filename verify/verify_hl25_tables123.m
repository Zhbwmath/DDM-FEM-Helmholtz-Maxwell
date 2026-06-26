% VERIFY_HL25_TABLES123  Reproduce Hu-Li Tables 1-3 under HPC gates.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'HL25_Helmholtz_harmonic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg = tableConfig();
cases = tableCases();
results = repmat(emptyResult(), numel(cases), 1);
fprintf('========== Hu-Li Tables 1-3 ==========\n');
fprintf('Cases: %d, permission gate: %.1f GB, runtime DOF cap: %d\n\n', ...
    numel(cases), cfg.permissionMemoryGB, cfg.maxRunDof);

if cfg.useParfor && cfg.runEnabled
    ensureParallelPool();
end

for i = 1:numel(cases)
    c = finalizeCase(cases(i));
    est = estimateCase(c, cfg);
    r = copyCase(emptyResult(), c, est);
    fprintf('[%02d/%02d] Table %s P%d eps=%s beta=%.1f k=%g*pi: N=%d, %.2f GB ... ', ...
        i, numel(cases), c.tableId, c.degree, c.epsilonLabel, ...
        c.beta, c.kappaPi, est.N, est.totalGB);

    if est.totalGB > cfg.hardMemoryGB
        r.status = 'blocked_memory_gt_hard_limit';
        r.notes = sprintf('Optimistic %s estimate %.2f GB exceeds hard limit %.2f GB.', ...
            est.memoryMode, est.totalGB, cfg.hardMemoryGB);
        fprintf('BLOCKED\n');
    elseif est.totalGB > cfg.permissionMemoryGB && ~cfg.allowPermissionRows
        r.status = 'requires_permission_gt_200gb';
        r.notes = sprintf('Optimistic %s estimate %.2f GB requires explicit permission.', ...
            est.memoryMode, est.totalGB);
        fprintf('PERMISSION REQUIRED\n');
    elseif ~cfg.runEnabled
        r.status = 'estimated_only';
        r.notes = 'Execution disabled by HL25_RUN=0.';
        fprintf('ESTIMATED\n');
    elseif ~cfg.runAllPermitted && est.N > cfg.maxRunDof
        r.status = 'queued_runtime_cap';
        r.notes = sprintf('Below memory gate but above runtime cap N=%d.', est.N);
        fprintf('QUEUED\n');
    else
        try
            runOut = runOneCase(c, cfg, est);
            r = mergeStruct(r, runOut);
            r.status = ternary(r.flag == 0, 'ran', 'ran_not_converged');
            r.iterationDifference = r.gmresIterations - r.paperIterations;
            r.ratioDifference = r.coarseRatio - r.paperRatio;
            fprintf('RAN (%d it, relres %.3e)\n', ...
                r.gmresIterations, r.finalRelres);
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
writeFigure(outDir, results);
fprintf('\nCSV: %s\n', fullfile(outDir, 'tables_1_3_results.csv'));
fprintf('Markdown: %s\n', fullfile(outDir, 'tables_1_3_results.md'));
fprintf('========== Hu-Li Tables 1-3 complete ==========\n');


function cfg = tableConfig()
cfg.permissionMemoryGB = envNumber('HL25_PERMISSION_MEMORY_GB', 200);
cfg.hardMemoryGB = envNumber('HL25_HARD_MEMORY_GB', 500);
cfg.allowPermissionRows = logical(envNumber('HL25_ALLOW_GT_200', 0));
cfg.maxRunDof = envNumber('HL25_MAX_RUN_DOF', 800000);
cfg.runAllPermitted = logical(envNumber('HL25_RUN_ALL_PERMITTED', 0));
cfg.runEnabled = logical(envNumber('HL25_RUN', 1));
cfg.useParfor = logical(envNumber('HL25_PARFOR', 0));
cfg.solverMode = envString('HL25_SOLVER_MODE', 'adaptive');
cfg.coarseSolverMode = envString('HL25_COARSE_SOLVER_MODE', 'adaptive');
cfg.rankMethod = envString('HL25_RANK_METHOD', 'none');
cfg.tol = envNumber('HL25_TOL', 1e-6);
cfg.maxit = envNumber('HL25_MAXIT', 1000);
cfg.restart = [];
cfg.localLuFillConstant = envNumber('HL25_LU_FILL', 20);
cfg.coarseLuFillConstant = envNumber('HL25_COARSE_LU_FILL', 20);
cfg.coarseNnzPerRow = envNumber('HL25_COARSE_NNZ_PER_ROW', 80);
cfg.gmresBasisLength = envNumber('HL25_GMRES_BASIS', 10);
cfg.timeLimitS = envNumber('HL25_TIME_LIMIT_S', 14400);
end


function cases = tableCases()
cases = repmat(baseCase(), 0, 1);

k1 = [30, 40, 50, 60];
cases = appendBlock(cases, '1', 1, 'k', [0.7, 0.6, 0.5], k1, ...
    [5 4 4; 4 3 3; 3 3 2; 3 2 2], ...
    [1.4 .5 .12; 2.8 .6 .14; 3.5 .9 .21; 4.5 1.0 .23]);
cases = appendBlock(cases, '1', 1, '0', [1.0, 0.8, 0.6], k1, ...
    [7 9 13; 6 6 10; 6 6 7; 4 5 6], ...
    [208.7 5.0 .23; 233.2 8.0 .24; 493.9 8.2 .28; 875.4 11.8 .29]);

k2 = [40, 60, 80, 100];
cases = appendBlock(cases, '1', 2, 'k', [0.7, 0.6, 0.5], k2, ...
    [6 4 3; 3 2 3; 3 2 2; 2 2 2], ...
    [3.7 1.1 .3; 11.9 2.1 .6; 18.7 3.7 .7; 25.9 4.4 1.0]);
cases = appendBlock(cases, '1', 2, '0', [1.0, 0.8, 0.6], k2, ...
    [6 6 10; 6 6 6; 4 4 5; 4 3 4], ...
    [507.8 22.4 .7; 1965.9 32.0 .9; 5064.8 73.4 1.4; 7714.3 113.6 1.5]);

k23 = [40, 80, 120, 160];
cases = appendBlock(cases, '2', 2, 'k', [0.7, 0.6, 0.5], k23, ...
    [7 6 6; 7 7 6; 8 7 7; 10 8 6], ...
    [3.3 .8 .15; 6.5 1.1 .18; 6.4 1.2 .17; 8.1 1.3 .17]);
cases = appendBlock(cases, '2', 2, '0', [0.7, 0.6, 0.5], k23, ...
    [7 6 6; 7 7 7; 9 8 8; 12 9 7], ...
    [3.3 .8 .15; 6.5 1.1 .18; 6.4 1.2 .17; 8.1 1.3 .17]);
cases = appendBlock(cases, '3', 2, 'k', [0.7, 0.6, 0.5], k23, ...
    [4 5 4; 5 4 5; 8 5 4; 7 6 5], ...
    [4.3 .9 .21; 7.2 1.4 .23; 7.2 1.4 .21; 9.8 1.5 .21]);
cases = appendBlock(cases, '3', 2, '0', [0.7, 0.6, 0.5], k23, ...
    [5 5 5; 5 5 5; 9 5 5; 7 6 5], ...
    [4.3 .9 .21; 7.2 1.4 .23; 7.2 1.4 .21; 9.8 1.5 .21]);
end


function cases = appendBlock(cases, tableId, degree, epsilonLabel, betas, ...
    kappaPi, iterations, ratios)
for i = 1:numel(kappaPi)
    for j = 1:numel(betas)
        c = baseCase();
        c.tableId = tableId;
        c.degree = degree;
        c.epsilonLabel = epsilonLabel;
        c.beta = betas(j);
        c.kappaPi = kappaPi(i);
        c.kappa = pi * kappaPi(i);
        c.paperIterations = iterations(i,j);
        c.paperRatio = ratios(i,j);
        if strcmp(tableId, '3')
            c.coarseType = 'economic';
            c.nu = max(1, round(c.kappa^(1 - c.beta)));
        else
            c.coarseType = 'spectral';
            if strcmp(tableId, '2')
                c.rho = 0.5 * c.kappa^((c.beta - 1) / 2);
            elseif strcmp(epsilonLabel, 'k')
                c.rho = c.kappa^(c.beta - 2 - 1 / (4 * degree));
            else
                c.rho = c.kappa^(-1);
            end
        end
        cases(end+1) = c; %#ok<AGROW>
    end
end
end


function c = baseCase()
c = struct('tableId', '', 'degree', NaN, 'epsilonLabel', '', ...
    'beta', NaN, 'kappaPi', NaN, 'kappa', NaN, 'coarseType', '', ...
    'rho', NaN, 'nu', NaN, 'paperIterations', NaN, 'paperRatio', NaN, ...
    'hInv', NaN, 'h', NaN, 'nSubSide', NaN, 'd', NaN, ...
    'deltaSteps', NaN, 'delta', NaN);
end


function c = finalizeCase(c)
rawSub = max(1, round(c.kappa^c.beta));
fineExponent = (2 * c.degree + 1) / (2 * c.degree);
rawFine = ceil(c.kappa^fineExponent);
c.nSubSide = rawSub;
c.hInv = ceil(rawFine / (4 * rawSub)) * (4 * rawSub);
c.h = 1 / c.hInv;
c.d = 1 / rawSub;
c.deltaSteps = c.hInv / (4 * rawSub);
c.delta = c.deltaSteps / c.hInv;
end


function est = estimateCase(c, cfg)
if c.degree == 1
    N = (c.hInv + 1)^2;
else
    N = (2 * c.hInv + 1)^2;
end
NT = 2 * c.hInv^2;
nSub = c.nSubSide^2;
localCells = ceil((c.d + 2 * c.delta) * c.hInv);
nl = (c.degree * localCells + 1)^2;
boundaryDof = max(4, 4 * c.degree * localCells);
if strcmp(c.coarseType, 'economic')
    rawCoarse = nSub * (2 * c.nu);
else
    rawCoarse = min(nSub * boundaryDof, ...
        max(1, ceil(c.paperRatio * nl)));
end
globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.gmresBasisLength + 3);
localMatrixBytes = nSub * 336 * nl;
storedLuBytes = nSub * 16 * cfg.localLuFillConstant * ...
    nl * log2(max(nl, 2));
directFactorPeakBytes = 16 * cfg.localLuFillConstant * ...
    nl * log2(max(nl, 2));
basisBytes = 16 * nl * rawCoarse;
[coarseBytes, coarseLuBytes, coarseDirectBytes, coarseSolverMode] = ...
    optimisticCoarseBytes(rawCoarse, cfg);
fixedBytes = globalBytes + gmresBytes + basisBytes + coarseBytes;
luBytes = fixedBytes + localMatrixBytes + storedLuBytes;
directBytes = fixedBytes + localMatrixBytes + directFactorPeakBytes;
[totalBytes, memoryMode, solverMode] = optimisticSolverBytes( ...
    luBytes, directBytes, cfg.solverMode, cfg.permissionMemoryGB);
est = struct('N', N, 'NT', NT, 'nSubdomains', nSub, ...
    'localDofEstimate', nl, 'boundaryDofEstimate', boundaryDof, ...
    'rawCoarseEstimate', rawCoarse, 'totalGB', totalBytes / 2^30, ...
    'luGB', luBytes / 2^30, 'directGB', directBytes / 2^30, ...
    'coarseLuGB', coarseLuBytes / 2^30, ...
    'coarseDirectGB', coarseDirectBytes / 2^30, ...
    'memoryMode', memoryMode, 'solverMode', solverMode, ...
    'coarseSolverMode', coarseSolverMode);
end


function out = runOneCase(c, cfg, est)
timer = tic;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
parts = partitionMesh2D(node, elem, bdFlag, ...
    [c.nSubSide, c.nSubSide], 'overlap', c.delta);
parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], ...
    [c.nSubSide, c.nSubSide], c.delta);
pde = helmholtzPDE(c.kappa, 'epsilon', epsilonValue(c.epsilonLabel));
opts = struct('degree', c.degree, 'coarseType', c.coarseType, ...
    'rho', c.rho, 'nu', c.nu, 'kappaRef', c.kappa, ...
    'solverMode', est.solverMode, 'useParfor', cfg.useParfor, ...
    'rankMethod', cfg.rankMethod, ...
    'coarseSolverMode', est.coarseSolverMode, ...
    'localLuFillConstant', cfg.localLuFillConstant);
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, parts, opts);
b = assemblePlaneWaveBoundaryLoad2D( ...
    node, elem, bdFlag, c.kappa, c.degree);
applyTime = 0;

    function y = timedApply(r)
        if toc(timer) > cfg.timeLimitS
            error('verify_hl25_tables123:timeLimit', ...
                'Case exceeded %.1f seconds.', cfg.timeLimitS);
        end
        applyTimer = tic;
        y = method.applyResidual(r);
        applyTime = applyTime + toc(applyTimer);
    end

[~, flag, relres, iter, resvec] = gmres( ...
    method.A, b, cfg.restart, cfg.tol, cfg.maxit, @timedApply);
out = struct('flag', flag, 'finalRelres', relres, ...
    'gmresIterations', gmresIterationCount(iter, cfg.restart), ...
    'resvecLength', numel(resvec), ...
    'coarseDimension', method.stats.coarseDimension, ...
    'rawCoarseDimension', method.stats.rawCoarseDimension, ...
    'coarseRatio', method.stats.coarseToMaxLocalRatio, ...
    'localDofMin', method.stats.localDofMin, ...
    'localDofMax', method.stats.localDofMax, ...
    'localDofMean', method.stats.localDofMean, ...
    'coarseSetupS', method.timing.totalSetup, ...
    'preconditionerApplyS', applyTime, 'totalS', toc(timer));
end


function r = emptyResult()
r = struct('tableId', '', 'degree', NaN, 'epsilonLabel', '', ...
    'beta', NaN, 'kappaPi', NaN, 'kappa', NaN, 'coarseType', '', ...
    'rho', NaN, 'nu', NaN, 'hInv', NaN, 'h', NaN, ...
    'nSubSide', NaN, 'nSubdomains', NaN, 'd', NaN, ...
    'deltaSteps', NaN, 'delta', NaN, 'ndof', NaN, ...
    'localDofEstimate', NaN, 'rawCoarseEstimate', NaN, ...
    'estimateGB', NaN, 'estimateLuGB', NaN, ...
    'estimateDirectGB', NaN, 'memoryMode', '', 'solverMode', '', ...
    'coarseSolverMode', '', ...
    'paperIterations', NaN, 'paperRatio', NaN, ...
    'gmresIterations', NaN, 'coarseDimension', NaN, ...
    'rawCoarseDimension', NaN, 'coarseRatio', NaN, ...
    'iterationDifference', NaN, 'ratioDifference', NaN, ...
    'flag', NaN, 'finalRelres', NaN, 'resvecLength', NaN, ...
    'localDofMin', NaN, 'localDofMax', NaN, 'localDofMean', NaN, ...
    'coarseSetupS', NaN, 'preconditionerApplyS', NaN, 'totalS', NaN, ...
    'status', '', 'notes', '');
end


function r = copyCase(r, c, est)
names = fieldnames(c);
for i = 1:numel(names)
    if isfield(r, names{i}), r.(names{i}) = c.(names{i}); end
end
r.ndof = est.N;
r.nSubdomains = est.nSubdomains;
r.localDofEstimate = est.localDofEstimate;
r.rawCoarseEstimate = est.rawCoarseEstimate;
r.estimateGB = est.totalGB;
r.estimateLuGB = est.luGB;
r.estimateDirectGB = est.directGB;
r.memoryMode = est.memoryMode;
r.solverMode = est.solverMode;
r.coarseSolverMode = est.coarseSolverMode;
end


function writeOutputs(outDir, results, cfg)
csvPath = fullfile(outDir, 'tables_1_3_results.csv');
mdPath = fullfile(outDir, 'tables_1_3_results.md');
writetable(struct2table(results), csvPath);

fid = fopen(mdPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Hu--Li, Tables 1--3.\n');
fprintf(fid, 'Created: 2026-06-10\n');
fprintf(fid, 'Updated: 2026-06-12\n');
fprintf(fid, 'Verification entry point: `verify/verify_hl25_tables123.m`\n');
fprintf(fid, 'Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `partitionMesh2D`, `linearPartitionOfUnity2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# Hu--Li Tables 1--3 Results\n\n');
fprintf(fid, 'Rows above %.1f GB require explicit permission; the default runtime cap is `N <= %d`. The constants in the article asymptotic relations are set to one before mesh alignment.\n\n', ...
    cfg.permissionMemoryGB, cfg.maxRunDof);
fprintf(fid, '| table | p | epsilon | beta | kappa | coarse | rho | nu | N | subdomains | estimate GB | memory mode | solver mode | paper it | repo it | paper ratio | repo ratio | status | notes |\n');
fprintf(fid, '|---|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %s | %d | %s | %.1f | %g pi | %s | %.3e | %g | %d | %d | %.2f | %s | %s | %g | %s | %g | %s | %s | %s |\n', ...
        r.tableId, r.degree, r.epsilonLabel, r.beta, r.kappaPi, ...
        r.coarseType, r.rho, r.nu, r.ndof, r.nSubdomains, ...
        r.estimateGB, r.memoryMode, r.solverMode, r.paperIterations, ...
        numberOrDash(r.gmresIterations), ...
        r.paperRatio, numberOrDash(r.coarseRatio), r.status, ...
        markdownEscape(r.notes));
end
end


function [bytes, luBytes, directBytes, solverMode] = optimisticCoarseBytes(n, cfg)
matrixBytes = 16 * max(n, cfg.coarseNnzPerRow * n);
luBytes = matrixBytes + 16 * cfg.coarseLuFillConstant * ...
    n * log2(max(n, 2));
directBytes = luBytes;
switch lower(cfg.coarseSolverMode)
    case {'lu', 'storedlu'}
        bytes = luBytes;
        solverMode = 'lu';
    case {'direct', 'backslash'}
        bytes = directBytes;
        solverMode = 'direct';
    case 'adaptive'
        if luBytes / 2^30 <= cfg.permissionMemoryGB
            bytes = luBytes;
            solverMode = 'lu';
        else
            bytes = directBytes;
            solverMode = 'direct';
        end
    otherwise
        error('verify_hl25_tables123:coarseSolverMode', ...
            'Unknown coarse solver mode "%s".', cfg.coarseSolverMode);
end
end


function [bytes, memoryMode, solverMode] = optimisticSolverBytes( ...
        luBytes, directBytes, requested, permissionGB)
switch lower(requested)
    case {'lu', 'storedlu'}
        bytes = luBytes;
        memoryMode = 'lu';
        solverMode = 'lu';
    case {'direct', 'backslash'}
        bytes = directBytes;
        memoryMode = 'direct';
        solverMode = 'direct';
    case 'adaptive'
        if directBytes < luBytes
            bytes = directBytes;
            memoryMode = 'direct';
        else
            bytes = luBytes;
            memoryMode = 'lu';
        end
        if luBytes / 2^30 <= permissionGB
            solverMode = 'lu';
        else
            solverMode = memoryMode;
        end
    otherwise
        error('verify_hl25_tables123:solverMode', ...
            'Unknown solver mode "%s".', requested);
end
end


function writeFigure(outDir, results)
ran = strcmp({results.status}, 'ran') | strcmp({results.status}, 'ran_not_converged');
if ~any(ran), return; end
r = results(ran);
figure('Name', 'Hu-Li Tables 1-3 comparison');
tiledlayout(1, 2);
nexttile;
plot([r.paperIterations], [r.gmresIterations], 'o', 'LineWidth', 1.2);
hold on;
limit = max([r.paperIterations, r.gmresIterations]);
plot([0, limit], [0, limit], 'k--', 'LineWidth', 1);
xlabel('Paper iterations', 'Interpreter', 'latex');
ylabel('Repository iterations', 'Interpreter', 'latex');
title('GMRES iteration comparison', 'Interpreter', 'latex');
grid on;
nexttile;
semilogy([r.paperRatio], [r.coarseRatio], 'o', 'LineWidth', 1.2);
hold on;
ratioLimit = max([r.paperRatio, r.coarseRatio]);
plot([min([r.paperRatio, r.coarseRatio]), ratioLimit], ...
    [min([r.paperRatio, r.coarseRatio]), ratioLimit], 'k--', 'LineWidth', 1);
xlabel('Paper coarse/local ratio', 'Interpreter', 'latex');
ylabel('Repository coarse/local ratio', 'Interpreter', 'latex');
title('Coarse-space size comparison', 'Interpreter', 'latex');
grid on;
exportgraphics(gcf, fullfile(outDir, 'tables_1_3_comparison.png'), 'Resolution', 180);
end


function out = mergeStruct(out, values)
names = fieldnames(values);
for i = 1:numel(names)
    out.(names{i}) = values.(names{i});
end
end


function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end


function n = gmresIterationCount(iter, restart)
if numel(iter) == 1
    n = iter;
elseif isempty(restart)
    n = iter(2);
else
    n = (iter(1) - 1) * restart + iter(2);
end
end


function text = numberOrDash(value)
if isnan(value), text = '-'; else, text = sprintf('%.6g', value); end
end


function text = markdownEscape(value)
text = strrep(char(value), '|', '\|');
text = strrep(text, newline, ' ');
end


function value = epsilonValue(label)
if strcmp(label, 'k'), value = 'k'; else, value = 0; end
end


function ensureParallelPool()
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', feature('numcores'));
end
end


function value = envNumber(name, defaultValue)
text = getenv(name);
if isempty(text)
    value = defaultValue;
else
    value = str2double(text);
    if isnan(value), value = defaultValue; end
end
end


function value = envString(name, defaultValue)
value = getenv(name);
if isempty(value), value = defaultValue; end
end
