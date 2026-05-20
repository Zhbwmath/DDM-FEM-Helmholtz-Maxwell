% VERIFY_ORAS_LARGEK_ITERATIONS  Large-k ORAS iteration tables.
%
%   Matrix-free Richardson/GMRES only.  No dense E is formed.
%   Resolution: h = 2*pi/(q*k), default q=10.
%   Strip overlap extension defaults to 1/4 on each side and can be
%   overridden with ORAS_LARGEK_STRIP_OVERLAP_EXTENSION.

fprintf('========== Large-k ORAS Iteration Study ==========\n\n');
opts = largeKOptions();
fprintf('Linear partition-of-unity weights are used for all overlap cases.\n');
fprintf('Memory budget for automatic runs: %.0f GB\n', opts.memoryBudgetGB);
fprintf('Strip overlap extension: %.4g (total adjacent overlap %.4g)\n', ...
    opts.stripOverlapExtension, 2*opts.stripOverlapExtension);
fprintf('Parpool mode: %s', opts.parpoolMode);
if ~isempty(opts.parpoolWorkers)
    fprintf(' (%d workers)', opts.parpoolWorkers);
end
fprintf('\n\n');

cases = buildCases(opts);
results = struct('shape', {}, 'k', {}, 'degree', {}, 'q', {}, 'hInv', {}, ...
    'dof', {}, 'grid', {}, 'memGB', {}, 'richardson', {}, 'gmres', {}, ...
    'relres', {}, 'flag', {}, 'seconds', {}, 'status', {});
outputDir = fileparts(mfilename('fullpath'));
matFile = taggedResultFile(fullfile(outputDir, 'oras_largek_iterations_results'), '.mat', opts.resultTag);
mdFile = taggedResultFile(fullfile(outputDir, 'oras_largek_iterations_results'), '.md', opts.resultTag);
fprintf('Checkpoint files: %s, %s\n', matFile, mdFile);
initMarkdownResults(mdFile, opts);

fprintf('%-10s %-5s %-3s %-4s %-7s %-10s %-8s %-10s %-12s %-8s %-8s\n', ...
    'shape', 'k', 'p', 'q', '1/h', 'DOF', 'grid', 'memGB', 'Richardson', 'GMRES', 'status');
fprintf('%s\n', repmat('-', 1, 105));

pool = maybeStartParpool(opts, numel(cases));
if isempty(pool)
    for c = 1:numel(cases)
        result = runConfiguredCase(cases(c), opts);
        results = checkpointResult(results, result, matFile, mdFile);
        printResult(result);
    end
else
    futures = parallel.FevalFuture.empty(0, numel(cases));
    for c = 1:numel(cases)
        futures(c) = parfeval(pool, @runConfiguredCase, 1, cases(c), opts);
    end
    for c = 1:numel(cases)
        [~, result] = fetchNext(futures);
        results = checkpointResult(results, result, matFile, mdFile);
        printResult(result);
    end
end

fprintf('\nResults written to %s\n', mdFile);
fprintf('========== Large-k ORAS iteration study complete ==========\n');


function opts = largeKOptions()
opts.kVals = numericEnvList('ORAS_LARGEK_KVALS', [40, 80, 120]);
opts.qVals = numericEnvList('ORAS_LARGEK_QVALS', 10);
opts.degrees = numericEnvList('ORAS_LARGEK_DEGREES', 1:3);
opts.parpoolMode = lower(strtrim(getenvDefault('ORAS_LARGEK_PARPOOL', 'off')));
opts.parpoolWorkers = numericEnvScalar('ORAS_LARGEK_WORKERS', []);
opts.memoryBudgetGB = numericEnvScalar('ORAS_LARGEK_MEMORY_GB', 400);
opts.stripOverlapExtension = numericEnvScalar('ORAS_LARGEK_STRIP_OVERLAP_EXTENSION', 1/4);
opts.resultTag = matlab.lang.makeValidName(strtrim(getenvDefault('ORAS_LARGEK_TAG', '')));
if strcmp(opts.resultTag, 'x')
    opts.resultTag = '';
end
if ~ismember(opts.parpoolMode, {'off', 'on', 'auto'})
    error('ORAS_LARGEK_PARPOOL must be off, on, or auto.');
end
end


function value = getenvDefault(name, defaultValue)
value = getenv(name);
if isempty(value)
    value = defaultValue;
end
end


function vals = numericEnvList(name, defaultValue)
txt = strtrim(getenv(name));
if isempty(txt)
    vals = defaultValue;
    return;
end
txt = strrep(txt, ',', ' ');
vals = str2num(txt); %#ok<ST2NM>
if isempty(vals)
    error('Environment variable %s must be a numeric list.', name);
end
end


function val = numericEnvScalar(name, defaultValue)
txt = strtrim(getenv(name));
if isempty(txt)
    val = defaultValue;
    return;
end
val = str2double(txt);
if isnan(val)
    error('Environment variable %s must be numeric.', name);
end
end


function fileName = taggedResultFile(stem, ext, tag)
if isempty(tag)
    fileName = [stem, ext];
else
    fileName = sprintf('%s_%s%s', stem, tag, ext);
end
end


function pool = maybeStartParpool(opts, nCases)
pool = [];
if strcmp(opts.parpoolMode, 'off') || nCases < 2
    return;
end
try
    if ~license('test', 'Distrib_Computing_Toolbox')
        if strcmp(opts.parpoolMode, 'on')
            error('Parallel Computing Toolbox license is unavailable.');
        end
        fprintf('Parallel Computing Toolbox unavailable; running serially.\n');
        return;
    end
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(opts.parpoolWorkers)
            pool = parpool('local');
        else
            pool = parpool('local', opts.parpoolWorkers);
        end
    end
catch err
    if strcmp(opts.parpoolMode, 'on')
        rethrow(err);
    end
    fprintf('Could not start parpool (%s); running serially.\n', err.message);
    pool = [];
end
end


function results = checkpointResult(results, result, matFile, mdFile)
results(end+1) = result; %#ok<AGROW>
save(matFile, 'results');
appendMarkdownResult(result, mdFile);
end


function result = runConfiguredCase(cfg, opts)
est = estimateCase(cfg);
result = struct('shape', cfg.shape, 'k', cfg.k, 'degree', cfg.degree, ...
    'q', cfg.q, 'hInv', round(1/cfg.h), 'dof', est.dof, ...
    'grid', cfg.gridLabel, 'memGB', est.totalGB, 'richardson', NaN, ...
    'gmres', NaN, 'relres', NaN, 'flag', NaN, 'seconds', NaN, 'status', "pending");

if est.totalGB > opts.memoryBudgetGB
    result.status = "skip-mem";
    return;
end

tStart = tic;
try
    [richIts, gmIts, relres, flag, dofActual] = runIterationCase(cfg);
    result.richardson = richIts;
    result.gmres = gmIts;
    result.relres = relres;
    result.flag = flag;
    result.dof = dofActual;
    result.seconds = toc(tStart);
    if flag == 0
        result.status = "ok";
    else
        result.status = "gmres-flag";
    end
catch err
    result.status = "failed";
    result.seconds = toc(tStart);
    fprintf('\nFAILED case %s k=%g p=%d: %s\n', cfg.shape, cfg.k, cfg.degree, err.message);
end
end


function cases = buildCases(opts)
kVals = opts.kVals;
qVals = opts.qVals;
degrees = opts.degrees;

cases = struct('shape', {}, 'k', {}, 'degree', {}, 'q', {}, 'h', {}, ...
    'nSub', {}, 'gridN', {}, 'gridLabel', {}, 'bbox', {}, 'overlap', {});

for k = kVals
    for p = degrees
        for q = qVals
            h = chooseResolutionMesh(k, q);
            cases(end+1) = stripCase(k, p, q, h, opts); %#ok<AGROW>
            cases(end+1) = gridCase(k, p, q, h); %#ok<AGROW>
        end
    end
end
end


function cfg = stripCase(k, degree, q, h, opts)
L = 16/3;
cfg = baseCase('strip', k, degree, q, h);
cfg.nSub = 8;
cfg.gridN = [];
cfg.gridLabel = '8x1';
cfg.bbox = [0, L, 0, 1];
cfg.overlap = opts.stripOverlapExtension;
cfg.h = alignedStripMeshSize(k, q, cfg.overlap);
end


function cfg = gridCase(k, degree, q, h)
gridN = max(2, round(k^0.4));
cfg = baseCase('grid', k, degree, q, h);
cfg.nSub = [];
cfg.gridN = gridN;
cfg.gridLabel = sprintf('%dx%d', gridN, gridN);
cfg.bbox = [0, 1, 0, 1];
cfg.overlap = 1/(4*gridN);
cfg.h = alignedGridMeshSize(k, q, gridN);
end


function cfg = baseCase(shape, k, degree, q, h)
cfg = struct();
cfg.shape = shape;
cfg.k = k;
cfg.degree = degree;
cfg.q = q;
cfg.h = h;
end


function h = chooseResolutionMesh(k, q)
hTarget = 2*pi/(q*k);
h = 1 / ceil(1/hTarget);
end


function h = alignedStripMeshSize(k, q, overlap)
% GGS strip tests use L=16/3 and H=2/3. For the common overlap extensions
% 1/4 and 1/2, choose h so both H/h and overlap/h are integers.
hTarget = 2*pi/(q*k);
if abs(overlap - 1/4) < 1e-12
    m = max(1, round(1/(12*hTarget)));
    h = 1/(12*m);
elseif abs(overlap - 1/2) < 1e-12
    m = max(1, round(1/(6*hTarget)));
    h = 1/(6*m);
else
    h = 1 / ceil(1/hTarget);
end
end


function h = alignedGridMeshSize(k, q, gridN)
% For checkerboards, H=1/gridN and overlap=H/4.  h=1/(4*gridN*m)
% resolves both the non-overlapping subdomains and the extended boundaries.
hTarget = 2*pi/(q*k);
m = max(1, round(1/(4*gridN*hTarget)));
h = 1/(4*gridN*m);
end


function budget = memoryBudgetGB()
budget = 400;
end


function est = estimateCase(cfg)
Lx = cfg.bbox(2) - cfg.bbox(1);
Ly = cfg.bbox(4) - cfg.bbox(3);
nx = ceil(Lx / cfg.h);
ny = ceil(Ly / cfg.h);

switch cfg.degree
    case 1
        dof = (nx+1) * (ny+1);
    case 2
        dof = (2*nx+1) * (2*ny+1);
    case 3
        dof = (3*nx+1) * (3*ny+1);
end

nnzPerRow = [9, 25, 49];
fillFactor = [12, 18, 25];
nnzA = nnzPerRow(cfg.degree) * dof;

if strcmp(cfg.shape, 'strip')
    nSub = cfg.nSub;
else
    nSub = cfg.gridN^2;
end

overlapFactor = estimateOverlapFactor(cfg);
localNnzTotal = overlapFactor * nnzA;

sparseGB = 32 * nnzA / 1e9;
luGB = 32 * fillFactor(cfg.degree) * localNnzTotal / 1e9;
krylovGB = 16 * dof * (maxGmresIter() + 10) / 1e9;
workGB = 4 * sparseGB + krylovGB;
totalGB = 1.5 * (sparseGB + luGB + workGB);

est = struct('dof', dof, 'nnzA', nnzA, 'nSub', nSub, ...
    'sparseGB', sparseGB, 'luGB', luGB, 'krylovGB', krylovGB, ...
    'totalGB', totalGB);
end


function fac = estimateOverlapFactor(cfg)
if strcmp(cfg.shape, 'strip')
    H = (cfg.bbox(2) - cfg.bbox(1)) / cfg.nSub;
    fac = 1 + 2*cfg.overlap/H;
else
    H = 1 / cfg.gridN;
    fac = (1 + 2*cfg.overlap/H)^2;
end
end


function maxIt = maxGmresIter()
maxIt = 120;
end


function [richIts, gmIts, relres, flag, dof] = runIterationCase(cfg)
[node, elem, bd, parts] = makeCase(cfg);
[nodeH, elemH] = extendForDegree(node, elem, cfg.degree);
[A, ~] = assembleHelmholtz2D(nodeH, elemH, bd, cfg.k, 0, 0, cfg.degree);
uExact = exp(1i * cfg.k * nodeH(:,1));
b = A * uExact;
dof = size(A, 1);

ap = orasHelmholtz(node, elem, bd, cfg.k, parts, cfg.degree);
richIts = richardsonIterations(A, b, ap, 1e-6, maxGmresIter());
[~, flag, relres, iter] = gmres(A, b, [], 1e-6, maxGmresIter(), ap);
gmIts = iter(2);
end


function [node, elem, bd, parts] = makeCase(cfg)
[node, elem, bd] = squaremesh(cfg.bbox, cfg.h);
if strcmp(cfg.shape, 'strip')
    parts = partitionMesh2D(node, elem, bd, cfg.nSub, 'overlap', cfg.overlap);
    parts = linearPartitionOfUnity2D(parts, cfg.bbox, [cfg.nSub, 1], cfg.overlap);
else
    parts = partitionMesh2D(node, elem, bd, [cfg.gridN, cfg.gridN], 'overlap', cfg.overlap);
    parts = linearPartitionOfUnity2D(parts, cfg.bbox, [cfg.gridN, cfg.gridN], cfg.overlap);
end
end


function [nodeH, elemH] = extendForDegree(node, elem, degree)
if degree == 1
    nodeH = node;
    elemH = elem;
else
    [nodeH, elemH] = extendMesh2D(node, elem, degree);
end
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


function printResult(r)
fprintf('%-10s %-5d %-3d %-4d %-7d %-10d %-8s %-10.1f %-12s %-8s %-8s\n', ...
    char(r.shape), r.k, r.degree, r.q, r.hInv, r.dof, char(r.grid), r.memGB, ...
    iterString(r.richardson), iterString(r.gmres), char(r.status));
end


function s = iterString(its)
if isnan(its)
    s = '-';
elseif isinf(its)
    s = 'skip';
elseif its > maxGmresIter()
    s = sprintf('>%d', maxGmresIter());
else
    s = sprintf('%d', its);
end
end


function initMarkdownResults(fileName, opts)
ensureParentDirectory(fileName);
fid = fopen(fileName, 'w');
if fid < 0
    error('Could not open %s for writing from %s.', fileName, pwd);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Large-k ORAS Iteration Results\n\n');
fprintf(fid, 'Resolution: `h = 2*pi/(q*k)`, currently `q=10`.\n\n');
fprintf(fid, 'Strip overlap extension: `%.8g` on each side; total adjacent overlap `%.8g`.\n\n', ...
    opts.stripOverlapExtension, 2*opts.stripOverlapExtension);
fprintf(fid, '| shape | k | p | q | 1/h | DOF | grid | est. GB | Richardson | GMRES | relres | flag | seconds | status |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n');
end


function appendMarkdownResult(r, fileName)
ensureParentDirectory(fileName);
fid = fopen(fileName, 'a');
if fid < 0
    error('Could not open %s for appending.', fileName);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '| %s | %d | %d | %d | %d | %d | %s | %.1f | %s | %s | %.2e | %.0f | %.1f | %s |\n', ...
    char(r.shape), r.k, r.degree, r.q, r.hInv, r.dof, char(r.grid), r.memGB, ...
    iterString(r.richardson), iterString(r.gmres), r.relres, r.flag, r.seconds, char(r.status));
end


function ensureParentDirectory(fileName)
[folder, ~, ~] = fileparts(fileName);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end
end
