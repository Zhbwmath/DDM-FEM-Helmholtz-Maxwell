% VERIFY_CIP_LXZZ_HULI_MEDIUM  CIP fine/local LXZZ with Hu-Li coarse spaces.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'Frmwrk_hybrd2lvl_DDM', 'CIP_Fine');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg = localConfig();
cases = localCases(cfg);
if cfg.resume
    results = loadResults(outDir, cfg);
else
    results = repmat(emptyRow(), 0, 1);
end

fprintf('========== CIP fine/local LXZZ with Hu-Li coarse space ==========\n\n');
fprintf('k values: %s\n', mat2str(cfg.kValues));
fprintf(['degree P%d, C_h=%g, beta=%g, Hu-Li overlap factor=%g, ' ...
    'run=%d, parfor=%d, workers=%d\n\n'], cfg.degree, cfg.Ch, ...
    cfg.beta, cfg.overlapFactor, cfg.runEnabled, cfg.useParfor, ...
    cfg.parpoolWorkers);

ensureParpool(cfg);

for i = 1:numel(cases)
    c = finalizeCase(cases(i), cfg);
    if cfg.resume && caseComplete(results, c)
        fprintf('[%02d/%02d] k=%g coarse=%s variant=%s: checkpoint complete, skipping.\n', ...
            i, numel(cases), c.k, c.coarseType, c.variant);
        continue;
    end

    est = estimateCase(c, cfg);
    row = baseRow(c, est, cfg);
    fprintf('[%02d/%02d] k=%g coarse=%s variant=%s: N=%d, coarse est=%d, %.2f GB ... ', ...
        i, numel(cases), c.k, c.coarseType, c.variant, est.N, ...
        est.rawCoarseEstimate, est.totalGB);

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
        row.notes = sprintf('Above run cap N=%d.', cfg.maxRunDof);
        fprintf('QUEUED DOF\n');
    elseif est.rawCoarseEstimate > cfg.maxCoarseEstimate
        row.status = 'queued_coarse_cap';
        row.notes = sprintf('Above coarse cap %d.', cfg.maxCoarseEstimate);
        fprintf('QUEUED COARSE\n');
    elseif ~cfg.runEnabled
        row.status = 'estimated_only';
        row.notes = 'Execution disabled by CIP_LXZZ_HULI_RUN=0.';
        fprintf('ESTIMATED\n');
    else
        row.status = 'running';
        row.notes = 'Checkpoint: row started.';
        results = replaceRow(results, row);
        writeOutputs(outDir, results, cfg);
        try
            row = runCase(row, c, est, cfg);
            fprintf('RAN: %d iterations, relres %.3e\n', ...
                row.gmresIterations, row.finalRelres);
        catch ME
            row.status = 'failed';
            row.notes = errorNote(ME);
            fprintf('FAILED: %s\n%s\n', ME.message, errorReport(ME));
        end
    end

    results = replaceRow(results, row);
    writeOutputs(outDir, results, cfg);
end

writeOutputs(outDir, results, cfg);
fprintf('\nCSV: %s\n', fullfile(outDir, [cfg.outputStem, '.csv']));
fprintf('Markdown: %s\n', fullfile(outDir, [cfg.outputStem, '.md']));
fprintf('========== CIP/Hu-Li LXZZ medium verification complete ==========\n');


function txt = errorNote(ME)
txt = sprintf('%s: %s', ME.identifier, ME.message);
if numel(txt) > 500
    txt = [txt(1:497), '...'];
end
end


function txt = errorReport(ME)
txt = getReport(ME, 'extended', 'hyperlinks', 'off');
end


function cfg = localConfig()
cfg.kValues = envVector('CIP_LXZZ_HULI_KVALUES', [16, 32, 64, 128]);
cfg.coarseTypes = envList('CIP_LXZZ_HULI_COARSE', {'economic'});
cfg.variants = envList('CIP_LXZZ_HULI_VARIANTS', {'dirichlet', 'impedance'});
cfg.degree = envNumber('CIP_LXZZ_HULI_DEGREE', 1);
cfg.beta = envNumber('CIP_LXZZ_HULI_BETA', 0.6);
cfg.Ch = envNumber('CIP_LXZZ_HULI_CH', 1);
cfg.gamma = envNumberOrEmpty('CIP_LXZZ_HULI_GAMMA', []);
cfg.cipOptions = struct();
cfg.overlapFactor = envNumber('CIP_LXZZ_HULI_OVERLAP_FACTOR', 1);
cfg.tol = envNumber('CIP_LXZZ_HULI_TOL', 1e-6);
cfg.maxit = envNumber('CIP_LXZZ_HULI_MAXIT', 100);
cfg.restart = [];
cfg.runEnabled = logical(envNumber('CIP_LXZZ_HULI_RUN', 1));
cfg.resume = logical(envNumber('CIP_LXZZ_HULI_RESUME', 1));
cfg.maxRunDof = envNumber('CIP_LXZZ_HULI_MAX_RUN_DOF', 50000);
cfg.maxCoarseEstimate = envNumber('CIP_LXZZ_HULI_MAX_COARSE', 300000);
cfg.permissionMemoryGB = envNumber('CIP_LXZZ_HULI_PERMISSION_GB', 200);
cfg.hardMemoryGB = envNumber('CIP_LXZZ_HULI_HARD_GB', 500);
cfg.allowPermissionRows = logical(envNumber('CIP_LXZZ_HULI_ALLOW_GT_200', 0));
cfg.solverMode = envString('CIP_LXZZ_HULI_SOLVER_MODE', 'lu');
cfg.huliSolverMode = envString('CIP_LXZZ_HULI_HULI_SOLVER_MODE', 'lu');
cfg.coarseSolverMode = envString('CIP_LXZZ_HULI_COARSE_SOLVER', 'lu');
cfg.adjointType = envString('CIP_LXZZ_HULI_ADJOINT', 'energy');
cfg.rankMethod = envString('CIP_LXZZ_HULI_RANK_METHOD', 'none');
cfg.useParfor = logical(envNumber('CIP_LXZZ_HULI_PARFOR', 1));
cfg.parpoolWorkers = envNumber('CIP_LXZZ_HULI_PARPOOL_WORKERS', feature('numcores'));
cfg.localSetupParfor = logical(envNumber('CIP_LXZZ_HULI_LOCAL_SETUP_PARFOR', 0));
cfg.localApplyParfor = logical(envNumber('CIP_LXZZ_HULI_LOCAL_APPLY_PARFOR', cfg.useParfor));
cfg.localApplyMode = envString('CIP_LXZZ_HULI_LOCAL_APPLY_MODE', 'auto');
cfg.fullVectorApplyLimitGB = envNumber('CIP_LXZZ_HULI_FULL_VECTOR_APPLY_LIMIT_GB', 2);
cfg.huliParfor = logical(envNumber('CIP_LXZZ_HULI_COARSE_PARFOR', cfg.useParfor));
cfg.localStorage = envString('CIP_LXZZ_HULI_LOCAL_STORAGE', 'matrix');
cfg.localLuFillConstant = envNumber('CIP_LXZZ_HULI_LU_FILL', 40);
cfg.huliLuFillConstant = envNumber('CIP_LXZZ_HULI_HULI_LU_FILL', 20);
cfg.coarseLuFillConstant = envNumber('CIP_LXZZ_HULI_COARSE_LU_FILL', 20);
cfg.gmresBasisLength = envNumber('CIP_LXZZ_HULI_GMRES_BASIS', 103);
cfg.outputStem = envString('CIP_LXZZ_HULI_OUTPUT_STEM', ...
    'cip_lxzz_huli_medium_results');
end


function ensureParpool(cfg)
needsPool = cfg.useParfor || cfg.localSetupParfor || ...
    cfg.localApplyParfor || cfg.huliParfor;
if ~needsPool
    return;
end
nWorkers = max(1, min(round(cfg.parpoolWorkers), feature('numcores')));
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting local parpool with %d workers ...\n', nWorkers);
    pool = parpool('local', nWorkers);
elseif pool.NumWorkers ~= nWorkers
    fprintf(['Existing parpool has %d workers; restarting with ' ...
        'requested %d workers ...\n'], pool.NumWorkers, nWorkers);
    delete(pool);
    pool = parpool('local', nWorkers);
else
    fprintf('Using existing parpool with requested %d workers.\n', ...
        pool.NumWorkers);
end
if pool.NumWorkers ~= nWorkers
    error('verify_cip_lxzz_huli_medium:parpoolWorkers', ...
        'Requested %d workers, active p.NumWorkers=%d.', ...
        nWorkers, pool.NumWorkers);
end
fprintf('Active parpool workers confirmed: p.NumWorkers=%d (requested %d).\n', ...
    pool.NumWorkers, nWorkers);
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
    'beta', NaN, 'rho', NaN, 'nu', NaN, 'hInv', NaN, ...
    'huliHInv', NaN, 'localHInv', NaN, 'h', NaN, 'huliH', NaN, ...
    'localH', NaN, 'overlap', NaN, 'huliNSubdomains', NaN, ...
    'localNSubdomains', NaN);
end


function c = finalizeCase(c, cfg)
c.degree = cfg.degree;
c.beta = cfg.beta;
c.rho = 0.5 * c.k^((c.beta - 1) / 2);
c.nu = max(1, round(c.k^(1 - c.beta)));
c.huliHInv = max(1, round(c.k));
rawFine = ceil(cfg.Ch * c.k^((2 * c.degree + 1) / (2 * c.degree)));
c.hInv = alignFineInv(rawFine, c.huliHInv);
c.h = 1 / c.hInv;
c.huliH = 1 / c.huliHInv;
overlapCells = max(1, round(cfg.overlapFactor * c.huliH * c.hInv));
c.overlap = overlapCells / c.hInv;
c.localH = selectLocalSpacing(c.k, c.variant);
c.localHInv = max(1, round(1 / c.localH));
c.huliNSubdomains = c.huliHInv^2;
c.localNSubdomains = (c.localHInv + 1)^2;
end


function spacing = selectLocalSpacing(k, variant)
switch lower(variant)
    case {'q1', 'dirichlet'}
        spacing = 1 / k;
    case {'q2', 'impedance'}
        spacing = 2 / k;
    otherwise
        error('verify_cip_lxzz_huli_medium:variant', ...
            'Unknown local solver variant "%s".', variant);
end
end


function est = estimateCase(c, cfg)
N = (c.degree * c.hInv + 1)^2;
NT = 2 * c.hInv^2;
lxzzSupportSteps = max(1, ceil(2 * c.hInv / c.localHInv));
lxzzLocalDof = (c.degree * lxzzSupportSteps + 1)^2;
huliSupportSteps = max(1, ceil((c.huliH + 2 * c.overlap) * c.hInv));
huliLocalDof = (c.degree * huliSupportSteps + 1)^2;
spectralBoundaryEstimate = max(2 * c.degree + 2, ...
    4 * c.degree * huliSupportSteps);
switch lower(c.coarseType)
    case 'economic'
        modesPerSubdomain = 2 * c.nu;
    case 'spectral'
        modesPerSubdomain = spectralBoundaryEstimate;
    otherwise
        error('verify_cip_lxzz_huli_medium:coarseType', ...
            'Unknown coarse type "%s".', c.coarseType);
end
rawCoarse = c.huliNSubdomains * modesPerSubdomain;

globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.gmresBasisLength + 3);
lxzzLocalMatrixBytes = c.localNSubdomains * 336 * lxzzLocalDof;
huliLocalMatrixBytes = c.huliNSubdomains * 336 * huliLocalDof;
huliLuBytes = c.huliNSubdomains * 16 * cfg.huliLuFillConstant * ...
    huliLocalDof * log2(max(huliLocalDof, 2));
huliBasisBytes = c.huliNSubdomains * 16 * huliLocalDof * modesPerSubdomain;
coarseBytes = rawCoarse * 80 * 16;
coarseLuBytes = 16 * cfg.coarseLuFillConstant * rawCoarse * ...
    log2(max(rawCoarse, 2));
energyBytes = 16 * 20 * N * log2(max(N, 2));
totalBytes = globalBytes + gmresBytes + lxzzLocalMatrixBytes + ...
    huliLocalMatrixBytes + huliLuBytes + huliBasisBytes + ...
    coarseBytes + coarseLuBytes + energyBytes;

est = struct('N', N, 'NT', NT, ...
    'lxzzSupportSteps', lxzzSupportSteps, ...
    'lxzzLocalDofEstimate', lxzzLocalDof, ...
    'huliSupportSteps', huliSupportSteps, ...
    'huliLocalDofEstimate', huliLocalDof, ...
    'modesPerSubdomain', modesPerSubdomain, ...
    'rawCoarseEstimate', rawCoarse, ...
    'totalGB', totalBytes / 2^30);
end


function row = runCase(row, c, est, cfg)
tAll = tic;
fprintf('\n    mesh setup ... ');
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
fprintf('done');

fprintf('\n    CIP fine space ... ');
fineOpts = struct('degree', c.degree, 'gamma', cfg.gamma, ...
    'cipOptions', cfg.cipOptions, ...
    'cacheEnergySolver', strcmpi(cfg.adjointType, 'energy'));
fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, c.k, fineOpts);
fprintf('done');

fprintf('\n    Hu-Li partition ... ');
huliParts = huliPartition(node, elem, c);
fprintf('done');

fprintf('\n    Hu-Li %s coarse basis ... ', c.coarseType);
pde = helmholtzPDE(c.k, 'epsilon', 0);
huliOpts = struct('degree', c.degree, 'coarseType', c.coarseType, ...
    'rho', c.rho, 'nu', c.nu, 'kappaRef', c.k, ...
    'solverMode', cfg.huliSolverMode, 'useParfor', cfg.huliParfor, ...
    'rankMethod', cfg.rankMethod, ...
    'coarseSolverMode', cfg.coarseSolverMode, ...
    'localLuFillConstant', cfg.huliLuFillConstant, ...
    'cacheEnergySolver', false, 'cacheEnergyAdjoint', false);
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, pde, huliParts, huliOpts);
huliStats = method.stats;
huliTiming = method.timing;
huliCoarse = method.coarseSpace;
clear huliParts
coarseEmbedding = embeddingForFineSpace(huliCoarse.trial, fine);
coarseSpace = struct('nativeTrial', huliCoarse.trial, ...
    'nativeTest', huliCoarse.test, 'embedding', coarseEmbedding, ...
    'description', sprintf(['Hu-Li %s Helmholtz-harmonic basis ' ...
    'with CIP fine Galerkin matrix'], lower(c.coarseType)), ...
    'coarseType', lower(c.coarseType), ...
    'rawDimension', huliCoarse.rawDimension, ...
    'dimension', huliCoarse.dimension, ...
    'rankInfo', huliCoarse.rankInfo);
clear method huliCoarse
fprintf('done (dim %d/%d)', huliStats.coarseDimension, ...
    huliStats.rawCoarseDimension);

fprintf('\n    CIP local partition ... ');
localParts = coarseHatPartition2D(node, elem, bdFlag, c.localH);
fprintf('done');

fprintf('\n    CIP local solver setup (%s storage) ... ', cfg.localStorage);
localOpts = struct('gamma', cfg.gamma, 'cipOptions', cfg.cipOptions, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'setupParfor', cfg.localSetupParfor, ...
    'applyParfor', cfg.localApplyParfor, ...
    'applyMode', cfg.localApplyMode, ...
    'fullVectorApplyLimitGB', cfg.fullVectorApplyLimitGB, ...
    'localLuFillConstant', cfg.localLuFillConstant, ...
    'localStorage', cfg.localStorage);
localSolver = buildCIPLxzzLocalSolversHelmholtz2D(fine, localParts, ...
    c.variant, localOpts);
fprintf('done');

fprintf('\n    hybrid assembly ... ');
preOpts = struct('fineSpace', fine, 'coarseSpace', coarseSpace, ...
    'localSolver', localSolver, 'variant', c.variant, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'adjointType', cfg.adjointType, ...
    'localLuFillConstant', cfg.localLuFillConstant);
pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, c.k, ...
    localParts, [], [], [], preOpts);
fprintf('done');

fprintf('\n    hybrid identity and GMRES ... ');
identity = verifyHybridIdentity(pre);
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.k, c.degree);
tSolve = tic;
[~, flag, relres, iter, resvec] = gmres( ...
    pre.A, b, cfg.restart, cfg.tol, cfg.maxit, @pre.applyResidual);

row.flag = flag;
row.finalRelres = relres;
row.gmresIterations = gmresIterationCount(iter, cfg.restart);
row.resvecLength = numel(resvec);
row.coarseDimension = size(pre.coarseSpace.AH, 1);
row.rawCoarseDimension = huliStats.rawCoarseDimension;
row.cipNonzeros = nnz(fine.cip.matrix);
row.solveS = toc(tSolve);
row.totalS = toc(tAll);
row.setupS = row.totalS - row.solveS;
row.huliSetupS = huliTiming.totalSetup;
row.localMode = pre.local.solverModeEffective;
if isfield(pre.local, 'localStorage') && ...
        ~strcmpi(pre.local.localStorage, 'factor')
    row.localMode = [pre.local.localStorage, '-', row.localMode];
end
if isfield(pre.local, 'applyMode')
    row.localMode = [row.localMode, '-', pre.local.applyMode];
end
row.localDofMin = pre.local.localDofMin;
row.localDofMax = pre.local.localDofMax;
row.localDofMean = pre.local.localDofMean;
row.identityError = identity;
row.status = ternary(flag == 0, 'ran', 'ran_not_converged');
row.notes = '';
if est.rawCoarseEstimate ~= huliStats.rawCoarseDimension
    row.notes = sprintf('Estimated raw coarse %d; actual %d.', ...
        est.rawCoarseEstimate, huliStats.rawCoarseDimension);
end
fprintf('done\n');
end


function err = verifyHybridIdentity(pre)
n = size(pre.A, 1);
rng(1);
x = randn(n, 1) + 1i * randn(n, 1);
y1 = pre.apply(x);
y2 = pre.applyResidual(pre.A * x);
err = norm(y1 - y2) / max(1, norm(y1));
if err > 1e-9
    error('verify_cip_lxzz_huli_medium:hybridIdentity', ...
        'apply(x) and applyResidual(A*x) differ by %.3e.', err);
end
end


function E = embeddingForFineSpace(nativeTrial, fine)
if size(nativeTrial, 1) == fine.N
    E = speye(fine.N);
elseif size(nativeTrial, 1) == size(fine.baseNode, 1)
    E = fine.p1ToFine;
else
    error('verify_cip_lxzz_huli_medium:coarseEmbedding', ...
        'Hu-Li coarse basis has %d rows but fine space has %d DOFs.', ...
        size(nativeTrial, 1), fine.N);
end
end


function parts = huliPartition(node, elem, c)
bbox = [0, 1, 0, 1];
gridSize = [c.huliHInv, c.huliHInv];
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


function row = emptyRow()
row = struct('k', NaN, 'degree', NaN, 'coarseType', '', ...
    'variant', '', 'beta', NaN, 'rho', NaN, 'nu', NaN, ...
    'fineHInv', NaN, 'huliHInv', NaN, 'localHInv', NaN, ...
    'overlap', NaN, 'ndof', NaN, 'huliSubdomains', NaN, ...
    'localSubdomains', NaN, 'huliSupportSteps', NaN, ...
    'localSupportSteps', NaN, 'huliLocalDofEstimate', NaN, ...
    'localDofEstimate', NaN, 'rawCoarseEstimate', NaN, ...
    'estimateGB', NaN, 'coarseDimension', NaN, ...
    'rawCoarseDimension', NaN, 'cipNonzeros', NaN, ...
    'localMode', '', 'localDofMin', NaN, 'localDofMax', NaN, ...
    'localDofMean', NaN, 'identityError', NaN, ...
    'gmresIterations', NaN, 'flag', NaN, 'finalRelres', NaN, ...
    'resvecLength', NaN, 'huliSetupS', NaN, 'setupS', NaN, ...
    'solveS', NaN, 'totalS', NaN, 'status', '', 'notes', '');
end


function row = baseRow(c, est, cfg)
row = emptyRow();
row.k = c.k;
row.degree = c.degree;
row.coarseType = c.coarseType;
row.variant = c.variant;
row.beta = c.beta;
row.rho = c.rho;
row.nu = c.nu;
row.fineHInv = c.hInv;
row.huliHInv = c.huliHInv;
row.localHInv = c.localHInv;
row.overlap = c.overlap;
row.ndof = est.N;
row.huliSubdomains = c.huliNSubdomains;
row.localSubdomains = c.localNSubdomains;
row.huliSupportSteps = est.huliSupportSteps;
row.localSupportSteps = est.lxzzSupportSteps;
row.huliLocalDofEstimate = est.huliLocalDofEstimate;
row.localDofEstimate = est.lxzzLocalDofEstimate;
row.rawCoarseEstimate = est.rawCoarseEstimate;
row.estimateGB = est.totalGB;
row.localMode = [cfg.localStorage, '-', cfg.solverMode];
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
fields = fieldnames(emptyRow());
for i = 1:height(t)
    r = emptyRow();
    for j = 1:numel(fields)
        name = fields{j};
        if ~ismember(name, t.Properties.VariableNames), continue; end
        v = t.(name)(i);
        if isnumeric(r.(name))
            r.(name) = double(v);
        elseif ismissing(v)
            r.(name) = '';
        else
            r.(name) = char(v);
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
if isempty(results), return; end
csvPath = fullfile(outDir, [cfg.outputStem, '.csv']);
mdPath = fullfile(outDir, [cfg.outputStem, '.md']);
writetable(struct2table(results, 'AsArray', true), csvPath);

fid = fopen(mdPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: CIP fine/local LXZZ with Hu--Li Helmholtz-harmonic coarse spaces.\n');
fprintf(fid, 'Created: 2026-07-02\n');
fprintf(fid, 'Updated: 2026-07-02\n');
fprintf(fid, 'Verification entry point: `verify/verify_cip_lxzz_huli_medium.m`\n');
fprintf(fid, 'Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# CIP Fine/Local LXZZ With Hu--Li Coarse Space\n\n');
fprintf(fid, 'Settings: P%d CIP fine and local forms, $h^{-1}=\\operatorname{align}(\\lceil C_h k^{(2p+1)/(2p)}\\rceil,k)$ with $C_h=%g$, Hu--Li $H=1/k$, Hu--Li overlap $\\delta=%gH$, $\\beta=%g$, local LXZZ spacing $1/k$ for Dirichlet and $2/k$ for impedance, adjoint `%s`, local storage `%s`, local apply mode `%s` with full-vector output limit %.2f GB, GMRES tolerance %.1e and max iterations %d.\n\n', ...
    cfg.degree, cfg.Ch, cfg.overlapFactor, cfg.beta, cfg.adjointType, ...
    cfg.localStorage, cfg.localApplyMode, cfg.fullVectorApplyLimitGB, ...
    cfg.tol, cfg.maxit);
fprintf(fid, '| k | coarse | variant | h^{-1} | H^{-1} | local H^{-1} | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |\n');
fprintf(fid, '|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %g | %s | %s | %g | %g | %g | %g | %s | %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
        r.k, r.coarseType, r.variant, r.fineHInv, r.huliHInv, ...
        r.localHInv, r.ndof, mdNumber(firstFinite(r.coarseDimension, ...
        r.rawCoarseEstimate)), r.estimateGB, mdNumber(r.cipNonzeros), ...
        r.localMode, mdNumber(r.localDofMax), mdNumber(r.identityError), ...
        mdNumber(r.gmresIterations), mdNumber(r.finalRelres), ...
        mdNumber(r.setupS), mdNumber(r.solveS), r.status, mdEscape(r.notes));
end
clear cleanup
end


function hInv = alignFineInv(raw, divisor)
hInv = max(1, ceil(raw));
if divisor > 0
    hInv = ceil(hInv / divisor) * divisor;
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


function n = envNumberOrEmpty(name, defaultValue)
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


function s = mdNumber(v)
if ~(isnumeric(v) && isscalar(v) && isfinite(v))
    s = '-';
elseif abs(v) >= 1e4 || (abs(v) > 0 && abs(v) < 1e-3)
    s = sprintf('%.3e', v);
else
    s = sprintf('%.6g', v);
end
end


function s = mdEscape(s)
s = char(string(s));
s = strrep(s, '|', '\|');
if isempty(s), s = ''; end
end
