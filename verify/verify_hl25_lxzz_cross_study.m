% VERIFY_HL25_LXZZ_CROSS_STUDY  Hu-Li coarse spaces in both LXZZ hybrids.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'HL25_Helmholtz_harmonic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg = crossConfig();
cases = crossCases(cfg);
results = loadCrossResults(outDir, cfg);
energySolverCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
cachedKappaPi = NaN;
fprintf('========== Hu-Li / LXZZ25 cross-study ==========\n\n');

if cfg.useParfor && cfg.runEnabled
    ensureParallelPool(cfg);
end

for i = 1:numel(cases)
    c = finalizeCrossCase(cases(i), cfg);
    if cfg.resume && caseAlreadyComplete(results, c)
        fprintf('[%02d/%02d] k=%g*pi eps=%s coarse=%s: checkpoint complete, skipping.\n', ...
            i, numel(cases), c.kappaPi, c.epsilonLabel, c.coarseType);
        continue;
    end
    if ~isequal(c.kappaPi, cachedKappaPi)
        energySolverCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
        cachedKappaPi = c.kappaPi;
    end
    est = estimateCrossCase(c, cfg);
    fprintf('[%02d/%02d] k=%g*pi eps=%s coarse=%s: N=%d, %.2f GB ... ', ...
        i, numel(cases), c.kappaPi, c.epsilonLabel, ...
        c.coarseType, est.N, est.totalGB);

    if est.totalGB > cfg.hardMemoryGB
        rows = classifiedRows(c, est, 'blocked_memory_gt_hard_limit', ...
            sprintf('Optimistic %s estimate %.2f GB exceeds hard limit.', ...
            est.memoryMode, est.totalGB));
        fprintf('BLOCKED\n');
    elseif est.totalGB > cfg.permissionMemoryGB && ~cfg.allowPermissionRows
        rows = classifiedRows(c, est, 'requires_permission_gt_200gb', ...
            sprintf('Optimistic %s estimate %.2f GB requires explicit permission.', ...
            est.memoryMode, est.totalGB));
        fprintf('PERMISSION REQUIRED\n');
    elseif ~cfg.runEnabled
        rows = classifiedRows(c, est, 'estimated_only', ...
            'Execution disabled by HL25_CROSS_RUN=0.');
        fprintf('ESTIMATED\n');
    elseif ~cfg.runAllPermitted && est.N > cfg.maxRunDof
        rows = classifiedRows(c, est, 'queued_runtime_cap', ...
            sprintf('Above relaxed runtime cap N=%d.', est.N));
        fprintf('QUEUED\n');
    else
        try
            rows = runCrossCase(c, est, cfg, energySolverCache, ...
                outDir, results);
            fprintf('RAN\n');
        catch ME
            rows = classifiedRows(c, est, 'failed', ...
                sprintf('%s: %s', ME.identifier, ME.message));
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    results = replaceCaseRows(results, c, rows);
    writeCrossOutputs(outDir, results, cfg);
end

writeCrossOutputs(outDir, results, cfg);
writeCrossFigure(outDir, results, cfg);
fprintf('\nCSV: %s\n', fullfile(outDir, [cfg.outputStem, '.csv']));
fprintf('Markdown: %s\n', fullfile(outDir, [cfg.outputStem, '.md']));
fprintf('========== Hu-Li / LXZZ25 cross-study complete ==========\n');


function cfg = crossConfig()
cfg.permissionMemoryGB = envNumber('HL25_CROSS_PERMISSION_MEMORY_GB', 200);
cfg.hardMemoryGB = envNumber('HL25_CROSS_HARD_MEMORY_GB', 500);
cfg.allowPermissionRows = logical(envNumber('HL25_CROSS_ALLOW_GT_200', 0));
cfg.resume = logical(envNumber('HL25_CROSS_RESUME', 0));
cfg.maxRunDof = envNumber('HL25_CROSS_MAX_RUN_DOF', 150000);
cfg.runAllPermitted = logical(envNumber('HL25_CROSS_RUN_ALL_PERMITTED', 0));
cfg.runEnabled = logical(envNumber('HL25_CROSS_RUN', 1));
cfg.useParfor = logical(envNumber('HL25_CROSS_PARFOR', 0));
cfg.parforWorkers = envNumber('HL25_CROSS_PARFOR_WORKERS', feature('numcores'));
cfg.solverMode = envString('HL25_CROSS_SOLVER_MODE', 'adaptive');
cfg.coarseSolverMode = envString('HL25_CROSS_COARSE_SOLVER_MODE', 'adaptive');
cfg.rankMethod = envString('HL25_CROSS_RANK_METHOD', 'none');
cfg.beta = envNumber('HL25_CROSS_BETA', 0.6);
cfg.lxzzLocalHFactor = envNumber('HL25_CROSS_LXZZ_H_FACTOR', 1);
if cfg.lxzzLocalHFactor <= 0
    error('verify_hl25_lxzz_cross_study:lxzzLocalHFactor', ...
        'HL25_CROSS_LXZZ_H_FACTOR must be positive.');
end
cfg.adjointType = envString('HL25_CROSS_ADJOINT', 'energy');
cfg.cacheEnergySolver = logical(envNumber('HL25_CROSS_CACHE_ENERGY_SOLVER', 1));
cfg.cacheEnergyAdjoint = logical(envNumber('HL25_CROSS_CACHE_ENERGY_ADJOINT', 0));
cfg.tol = envNumber('HL25_CROSS_TOL', 1e-6);
cfg.maxit = envNumber('HL25_CROSS_MAXIT', 1000);
cfg.restart = [];
cfg.localLuFillConstant = envNumber('HL25_CROSS_LU_FILL', 20);
cfg.coarseLuFillConstant = envNumber('HL25_CROSS_COARSE_LU_FILL', 20);
cfg.coarseNnzPerRow = envNumber('HL25_CROSS_COARSE_NNZ_PER_ROW', 80);
cfg.gmresBasisLength = envNumber('HL25_CROSS_GMRES_BASIS', 10);
cfg.spectralCoarseRatio = envNumber('HL25_CROSS_SPECTRAL_RATIO', 1.0);
cfg.energyFillConstant = envNumber('HL25_CROSS_ENERGY_FILL', 10);
cfg.timeLimitS = envNumber('HL25_CROSS_TIME_LIMIT_S', 14400);
cfg.outputStem = envString('HL25_CROSS_OUTPUT_STEM', 'lxzz_cross_results');
cfg.figureStem = envString('HL25_CROSS_FIGURE_STEM', 'lxzz_cross_iterations');
cfg.kappaPi = envNumberList('HL25_CROSS_KAPPA_PI', ...
    [8, 16, 40, 80, 120, 160]);
end


function cases = crossCases(cfg)
kappaPi = cfg.kappaPi;
epsilonLabels = {'0', 'k'};
coarseTypes = {'economic', 'spectral'};
cases = repmat(baseCrossCase(), 0, 1);
for i = 1:numel(kappaPi)
    for j = 1:numel(epsilonLabels)
        for q = 1:numel(coarseTypes)
            c = baseCrossCase();
            c.kappaPi = kappaPi(i);
            c.kappa = pi * kappaPi(i);
            c.epsilonLabel = epsilonLabels{j};
            c.coarseType = coarseTypes{q};
            c.beta = cfg.beta;
            c.rho = 0.5 * c.kappa^((c.beta - 1) / 2);
            c.nu = max(1, round(c.kappa^(1 - c.beta)));
            cases(end+1) = c; %#ok<AGROW>
        end
    end
end
end


function c = baseCrossCase()
c = struct('kappaPi', NaN, 'kappa', NaN, 'epsilonLabel', '', ...
    'coarseType', '', 'beta', NaN, 'rho', NaN, 'nu', NaN, ...
    'degree', 2, 'hInv', NaN, 'h', NaN, 'nSubSide', NaN, ...
    'd', NaN, 'deltaSteps', NaN, 'delta', NaN, ...
    'lxzzNSubSide', NaN, 'lxzzNSubdomains', NaN, ...
    'lxzzH', NaN, 'lxzzDeltaSteps', NaN, 'lxzzDelta', NaN);
end


function c = finalizeCrossCase(c, cfg)
c.nSubSide = max(1, round(c.kappa^c.beta));
rawFine = ceil(c.kappa^((2 * c.degree + 1) / (2 * c.degree)));
c.hInv = ceil(rawFine / (4 * c.nSubSide)) * (4 * c.nSubSide);
c.h = 1 / c.hInv;
c.d = 1 / c.nSubSide;
c.deltaSteps = c.hInv / (4 * c.nSubSide);
c.delta = c.deltaSteps / c.hInv;
lxzzTarget = max(1, ceil(c.kappa / cfg.lxzzLocalHFactor));
c.lxzzNSubSide = alignedSubdomainDivisor(c.hInv, lxzzTarget);
c.lxzzNSubdomains = (c.lxzzNSubSide + 1)^2;
c.lxzzH = 1 / c.lxzzNSubSide;
c.lxzzDeltaSteps = max(1, round(c.hInv / c.lxzzNSubSide));
c.lxzzDelta = c.lxzzDeltaSteps / c.hInv;
end


function nSubSide = alignedSubdomainDivisor(hInv, target)
if target <= 1
    nSubSide = 1;
    return;
end
divisors = find(mod(hInv, 1:hInv) == 0);
candidates = divisors(divisors >= target);
if isempty(candidates)
    nSubSide = hInv;
else
    nSubSide = candidates(1);
end
end


function est = estimateCrossCase(c, cfg)
N = (2 * c.hInv + 1)^2;
NT = 2 * c.hInv^2;
nSub = c.nSubSide^2;
localCells = ceil((c.d + 2 * c.delta) * c.hInv);
nl = (2 * localCells + 1)^2;
lxzzSupportSteps = max(1, ceil(2 * c.hInv / c.lxzzNSubSide));
lxzzNl = (2 * lxzzSupportSteps + 1)^2;
boundaryDof = max(4, 8 * localCells);
if strcmp(c.coarseType, 'economic')
    rawCoarse = nSub * (2 * c.nu);
else
    rawCoarse = min(nSub * boundaryDof, ...
        max(1, ceil(cfg.spectralCoarseRatio * nl)));
end
globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.gmresBasisLength + 3);
localMatrixBytes = nSub * 336 * nl;
storedLuBytes = nSub * 16 * cfg.localLuFillConstant * ...
    nl * log2(max(nl, 2));
directFactorPeakBytes = 16 * cfg.localLuFillConstant * ...
    nl * log2(max(nl, 2));
lxzzLocalMatrixBytes = c.lxzzNSubdomains * 336 * lxzzNl;
lxzzStoredLuBytes = c.lxzzNSubdomains * 16 * cfg.localLuFillConstant * ...
    lxzzNl * log2(max(lxzzNl, 2));
lxzzDirectFactorPeakBytes = 16 * cfg.localLuFillConstant * ...
    lxzzNl * log2(max(lxzzNl, 2));
builderBasisBytes = 16 * nl * rawCoarse;
[coarseBytes, coarseLuBytes, coarseDirectBytes, coarseSolverMode] = ...
    optimisticCoarseBytes(rawCoarse, cfg);
energyFactorBytes = 0;
if strcmpi(cfg.adjointType, 'energy')
    energyFactorBytes = 16 * cfg.energyFillConstant * N * log2(max(N, 2));
end
energyAdjointBytes = 0;
if cfg.cacheEnergyAdjoint
    energyAdjointBytes = 16 * N * rawCoarse;
end
fixedBytes = globalBytes + gmresBytes + builderBasisBytes + ...
    coarseBytes + energyFactorBytes + energyAdjointBytes;
luBytes = fixedBytes + localMatrixBytes + storedLuBytes + ...
    lxzzLocalMatrixBytes + lxzzStoredLuBytes;
directBytes = fixedBytes + localMatrixBytes + lxzzLocalMatrixBytes + ...
    directFactorPeakBytes + lxzzDirectFactorPeakBytes;
[totalBytes, memoryMode, solverMode] = optimisticSolverBytes( ...
    luBytes, directBytes, cfg.solverMode, cfg.permissionMemoryGB);
est = struct('N', N, 'NT', NT, 'nSubdomains', nSub, ...
    'localDofEstimate', nl, 'rawCoarseEstimate', rawCoarse, ...
    'lxzzNSubdomains', c.lxzzNSubdomains, ...
    'lxzzLocalDofEstimate', lxzzNl, ...
    'totalGB', totalBytes / 2^30, 'luGB', luBytes / 2^30, ...
    'directGB', directBytes / 2^30, ...
    'coarseLuGB', coarseLuBytes / 2^30, ...
    'coarseDirectGB', coarseDirectBytes / 2^30, ...
    'memoryMode', memoryMode, 'solverMode', solverMode, ...
    'coarseSolverMode', coarseSolverMode);
end


function rows = runCrossCase(c, est, cfg, energySolverCache, outDir, priorResults)
timer = tic;
fprintf('\n    assembling mesh, Hu-Li partition, LXZZ local partition, and Hu-Li coarse space ...\n');
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
hlParts = partitionMesh2D(node, elem, bdFlag, ...
    [c.nSubSide, c.nSubSide], 'overlap', c.delta);
hlParts = linearPartitionOfUnity2D(hlParts, [0, 1, 0, 1], ...
    [c.nSubSide, c.nSubSide], c.delta);
lxzzParts = coarseHatPartition2D(node, elem, bdFlag, 1 / c.lxzzNSubSide);
pde = helmholtzPDE(c.kappa, 'epsilon', epsilonValue(c.epsilonLabel));
builderOpts = struct('degree', 2, 'coarseType', c.coarseType, ...
    'rho', c.rho, 'nu', c.nu, 'kappaRef', c.kappa, ...
    'solverMode', est.solverMode, 'useParfor', cfg.useParfor, ...
    'rankMethod', cfg.rankMethod, ...
    'coarseSolverMode', est.coarseSolverMode, ...
    'localLuFillConstant', cfg.localLuFillConstant, ...
    'cacheEnergySolver', cfg.cacheEnergySolver && ...
        strcmpi(cfg.adjointType, 'energy'), ...
    'cacheEnergyAdjoint', cfg.cacheEnergyAdjoint);
energyKey = sprintf('k%.16g_p%d_h%d', c.kappa, c.degree, c.hInv);
if isKey(energySolverCache, energyKey)
    builderOpts.energySolve = energySolverCache(energyKey);
end
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, hlParts, builderOpts);
fprintf('    Hu-Li setup %.2f s, coarse dim %d/%d ...\n', ...
    method.timing.totalSetup, method.stats.coarseDimension, ...
    method.stats.rawCoarseDimension);
if cfg.cacheEnergySolver && strcmpi(cfg.adjointType, 'energy') && ...
        ~isKey(energySolverCache, energyKey)
    energySolverCache(energyKey) = method.fineSpace.energySolve; %#ok<NASGU>
end
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.kappa, 2);

rows = repmat(baseRunRow(c, est, method), 3, 1);
rows(1).method = 'hu_li_native';
rows(2).method = 'lxzz_dirichlet';
rows(3).method = 'lxzz_impedance';
for j = 1:3
    rows(j).status = 'running_pending';
    rows(j).notes = 'Group started; method has not completed yet.';
end
checkpointMethodRows(outDir, priorResults, c, rows, cfg);
fprintf('    solving Hu-Li native ...\n');
rows(1) = solveOne(rows(1), 'hu_li_native', method.A, b, ...
    method.applyResidual, method.timing.totalSetup, timer, cfg);
checkpointMethodRows(outDir, priorResults, c, rows, cfg);
fprintf('    Hu-Li native done: %d iterations, relres %.3e ...\n', ...
    rows(1).gmresIterations, rows(1).finalRelres);

variants = {'dirichlet', 'impedance'};
names = {'lxzz_dirichlet', 'lxzz_impedance'};
for j = 1:2
    fprintf('    building and solving %s ...\n', names{j});
    preOpts = struct('fineSpace', method.fineSpace, ...
        'coarseSpace', method.coarseSpace, 'variant', variants{j}, ...
        'solverMode', est.solverMode, 'useParfor', cfg.useParfor, ...
        'adjointType', cfg.adjointType, ...
        'localLuFillConstant', cfg.localLuFillConstant);
    precon = twoLevelHybridSchwarzHelmholtz2D( ...
        node, elem, bdFlag, pde, lxzzParts, [], [], [], preOpts);
    rows(j+1) = solveOne(rows(j+1), names{j}, precon.A, b, ...
        precon.applyResidual, method.timing.totalSetup + ...
        precon.timing.localSetup, timer, cfg);
    checkpointMethodRows(outDir, priorResults, c, rows, cfg);
    fprintf('    %s done: %d iterations, relres %.3e ...\n', ...
        names{j}, rows(j+1).gmresIterations, rows(j+1).finalRelres);
end
end


function checkpointMethodRows(outDir, priorResults, c, rows, cfg)
checkpoint = replaceCaseRows(priorResults, c, rows);
writeCrossOutputs(outDir, checkpoint, cfg);
end


function row = solveOne(row, methodName, A, b, applyInverse, setupTime, timer, cfg)
applyTime = 0;
    function y = timedApply(r)
        if toc(timer) > cfg.timeLimitS
            error('verify_hl25_lxzz_cross_study:timeLimit', ...
                'Cross-study group exceeded %.1f seconds.', cfg.timeLimitS);
        end
        applyTimer = tic;
        y = applyInverse(r);
        applyTime = applyTime + toc(applyTimer);
    end

solveTimer = tic;
[~, flag, relres, iter, resvec] = gmres( ...
    A, b, cfg.restart, cfg.tol, cfg.maxit, @timedApply);
row.method = methodName;
row.flag = flag;
row.finalRelres = relres;
row.gmresIterations = gmresIterationCount(iter, cfg.restart);
row.resvecLength = numel(resvec);
row.setupS = setupTime;
row.applyS = applyTime;
row.solveS = toc(solveTimer);
row.totalGroupS = toc(timer);
row.status = ternary(flag == 0, 'ran', 'ran_not_converged');
row.notes = '';
end


function rows = classifiedRows(c, est, status, notes)
names = {'hu_li_native', 'lxzz_dirichlet', 'lxzz_impedance'};
rows = repmat(emptyCrossResult(), 3, 1);
for j = 1:3
    rows(j) = baseRunRow(c, est, []);
    rows(j).method = names{j};
    rows(j).status = status;
    rows(j).notes = notes;
end
end


function row = baseRunRow(c, est, method)
row = emptyCrossResult();
names = fieldnames(c);
for i = 1:numel(names)
    if isfield(row, names{i}), row.(names{i}) = c.(names{i}); end
end
row.ndof = est.N;
row.nSubdomains = est.nSubdomains;
row.localDofEstimate = est.localDofEstimate;
row.rawCoarseEstimate = est.rawCoarseEstimate;
row.lxzzNSubdomains = est.lxzzNSubdomains;
row.lxzzLocalDofEstimate = est.lxzzLocalDofEstimate;
row.estimateGB = est.totalGB;
row.memoryMode = est.memoryMode;
row.solverMode = est.solverMode;
row.coarseSolverMode = est.coarseSolverMode;
row.estimateLuGB = est.luGB;
row.estimateDirectGB = est.directGB;
if ~isempty(method)
    row.coarseDimension = method.stats.coarseDimension;
    row.rawCoarseDimension = method.stats.rawCoarseDimension;
    row.coarseRatio = method.stats.coarseToMaxLocalRatio;
end
end


function row = emptyCrossResult()
row = struct('method', '', 'kappaPi', NaN, 'kappa', NaN, ...
    'epsilonLabel', '', 'coarseType', '', 'beta', NaN, ...
    'rho', NaN, 'nu', NaN, 'degree', NaN, 'hInv', NaN, 'h', NaN, ...
    'nSubSide', NaN, 'nSubdomains', NaN, 'd', NaN, ...
    'deltaSteps', NaN, 'delta', NaN, ...
    'lxzzNSubSide', NaN, 'lxzzNSubdomains', NaN, ...
    'lxzzH', NaN, 'lxzzDeltaSteps', NaN, 'lxzzDelta', NaN, ...
    'ndof', NaN, ...
    'localDofEstimate', NaN, 'rawCoarseEstimate', NaN, ...
    'lxzzLocalDofEstimate', NaN, ...
    'estimateGB', NaN, 'estimateLuGB', NaN, ...
    'estimateDirectGB', NaN, 'memoryMode', '', 'solverMode', '', ...
    'coarseSolverMode', '', ...
    'coarseDimension', NaN, ...
    'rawCoarseDimension', NaN, 'coarseRatio', NaN, ...
    'gmresIterations', NaN, 'flag', NaN, 'finalRelres', NaN, ...
    'resvecLength', NaN, 'setupS', NaN, 'applyS', NaN, ...
    'solveS', NaN, 'totalGroupS', NaN, 'status', '', 'notes', '');
end


function writeCrossOutputs(outDir, results, cfg)
csvPath = fullfile(outDir, [cfg.outputStem, '.csv']);
mdPath = fullfile(outDir, [cfg.outputStem, '.md']);
writetable(struct2table(results), csvPath);

fid = fopen(mdPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.\n');
fprintf(fid, 'Created: 2026-06-10\n');
fprintf(fid, 'Updated: 2026-06-22\n');
fprintf(fid, 'Verification entry point: `verify/verify_hl25_lxzz_cross_study.m`\n');
fprintf(fid, 'Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# Hu--Li / LXZZ25 Cross-Study Results\n\n');
fprintf(fid, 'Every coarse-space configuration is tested with the native Hu--Li hybrid, LXZZ Dirichlet hybrid, and LXZZ impedance hybrid for both $\\epsilon=0$ and $\\epsilon=\\kappa$. Hu--Li coarse spaces use the Hu--Li paper partition; LXZZ local solvers use a separate aligned `coarseHatPartition2D` partition with $H_{\\rm LXZZ}\\le 1/\\kappa$ by default. Rows above %.1f GB require explicit permission.\n\n', ...
    cfg.permissionMemoryGB);
fprintf(fid, '| method | coarse | epsilon | kappa | N | Hu-Li subdomains | LXZZ local subdomains | coarse dim | ratio | estimate GB | memory mode | local solver | coarse solver | iterations | relres | setup s | solve s | status | notes |\n');
fprintf(fid, '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %s | %s | %s | %g pi | %d | %d | %d | %s | %s | %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
        r.method, r.coarseType, r.epsilonLabel, r.kappaPi, r.ndof, ...
        r.nSubdomains, r.lxzzNSubdomains, numberOrDash(r.coarseDimension), ...
        numberOrDash(r.coarseRatio), r.estimateGB, r.memoryMode, ...
        r.solverMode, r.coarseSolverMode, ...
        numberOrDash(r.gmresIterations), numberOrDash(r.finalRelres), ...
        numberOrDash(r.setupS), numberOrDash(r.solveS), ...
        r.status, markdownEscape(r.notes));
end
end


function results = loadCrossResults(outDir, cfg)
results = repmat(emptyCrossResult(), 0, 1);
if ~cfg.resume
    return;
end
csvPath = fullfile(outDir, [cfg.outputStem, '.csv']);
if ~exist(csvPath, 'file')
    return;
end
importOpts = detectImportOptions(csvPath, 'TextType', 'string');
textVars = {'method', 'epsilonLabel', 'coarseType', 'memoryMode', ...
    'solverMode', 'coarseSolverMode', 'status', 'notes'};
presentTextVars = intersect(textVars, importOpts.VariableNames);
if ~isempty(presentTextVars)
    importOpts = setvartype(importOpts, presentTextVars, 'string');
end
T = readtable(csvPath, importOpts);
raw = table2struct(T);
for i = 1:numel(raw)
    row = emptyCrossResult();
    names = fieldnames(row);
    for j = 1:numel(names)
        if isfield(raw, names{j})
            value = raw(i).(names{j});
            if isstring(value), value = char(value); end
            if strcmp(names{j}, 'epsilonLabel') && isnumeric(value)
                value = num2str(value);
            end
            row.(names{j}) = value;
        end
    end
    results(end+1) = row; %#ok<AGROW>
end
end


function yes = caseAlreadyComplete(results, c)
mask = crossCaseMask(results, c);
if nnz(mask) < 3
    yes = false;
    return;
end
status = {results(mask).status};
terminal = {'ran', 'ran_not_converged', ...
    'blocked_memory_gt_hard_limit', 'requires_permission_gt_200gb', ...
    'queued_runtime_cap'};
yes = all(ismember(status, terminal));
end


function results = replaceCaseRows(results, c, rows)
mask = crossCaseMask(results, c);
old = results(~mask);
results = repmat(emptyCrossResult(), numel(old) + numel(rows), 1);
for i = 1:numel(old)
    results(i) = old(i);
end
for j = 1:numel(rows)
    results(numel(old) + j) = rows(j);
end
end


function mask = crossCaseMask(results, c)
if isempty(results)
    mask = false(0, 1);
    return;
end
mask = abs([results.kappaPi] - c.kappaPi) < 1e-12 & ...
    abs([results.beta] - c.beta) < 1e-12 & ...
    strcmp({results.epsilonLabel}, c.epsilonLabel) & ...
    strcmp({results.coarseType}, c.coarseType);
mask = mask(:);
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
        error('verify_hl25_lxzz_cross_study:coarseSolverMode', ...
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
        error('verify_hl25_lxzz_cross_study:solverMode', ...
            'Unknown solver mode "%s".', requested);
end
end


function writeCrossFigure(outDir, results, cfg)
ran = strcmp({results.status}, 'ran') | strcmp({results.status}, 'ran_not_converged');
if ~any(ran), return; end
r = results(ran);
methods = unique({r.method}, 'stable');
figure('Name', 'Hu-Li / LXZZ25 cross-study');
hold on;
for j = 1:numel(methods)
    use = strcmp({r.method}, methods{j});
    plot([r(use).kappaPi], [r(use).gmresIterations], '-o', ...
        'LineWidth', 1.2, 'DisplayName', strrep(methods{j}, '_', '\_'));
end
xlabel('$\kappa/\pi$', 'Interpreter', 'latex');
ylabel('GMRES iterations', 'Interpreter', 'latex');
title('Hu--Li coarse spaces in hybrid Schwarz methods', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;
exportgraphics(gcf, fullfile(outDir, [cfg.figureStem, '.png']), 'Resolution', 180);
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


function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
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


function ensureParallelPool(cfg)
pool = gcp('nocreate');
if isempty(pool)
    workers = max(1, min(feature('numcores'), round(cfg.parforWorkers)));
    parpool('local', workers);
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


function values = envNumberList(name, defaultValues)
text = getenv(name);
if isempty(text)
    values = defaultValues;
    return;
end
pieces = regexp(text, '[,; ]+', 'split');
values = [];
for i = 1:numel(pieces)
    if isempty(pieces{i}), continue; end
    v = str2double(pieces{i});
    if ~isnan(v)
        values(end+1) = v; %#ok<AGROW>
    end
end
if isempty(values)
    values = defaultValues;
end
end
