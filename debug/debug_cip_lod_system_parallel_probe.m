% DEBUG_CIP_LOD_SYSTEM_PARALLEL_PROBE  Build only the LOD system with parfor.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

outDir = fullfile(repoRoot, 'debug');
if ~exist(outDir, 'dir'), mkdir(outDir); end

k = envNumber('CIP_LOD_PROBE_K', 128);
ell = envNumber('CIP_LOD_PROBE_ELL', 2);
workers = envNumber('CIP_LOD_PROBE_WORKERS', 4);
Ch = envNumber('CIP_LOD_PROBE_CH', 1);
degree = 1;
HInv = max(1, round(k));
rawFine = ceil(Ch * k^((2 * degree + 1) / (2 * degree)));
hInv = ceil(rawFine / HInv) * HInv;
h = 1 / hInv;
H = 1 / HInv;

stamp = datestr(now, 'yyyymmdd_HHMMSS');
outBase = sprintf('cip_lod_system_probe_k%d_hinv%d_Hinv%d_ell%d_w%d_%s', ...
    round(k), hInv, HInv, ell, workers, stamp);
matPath = fullfile(outDir, [outBase, '.mat']);
logPath = fullfile(outDir, [outBase, '.log']);
diary(logPath);
cleanupDiary = onCleanup(@() diary('off'));

fprintf('========== CIP LOD system parallel probe ==========\n');
fprintf('k=%g, h^{-1}=%d, H^{-1}=%d, ell=%d, requested workers=%d\n', ...
    k, hInv, HInv, ell, workers);
fprintf('Started: %s\n', datestr(now, 31));

checkpoint('config', matPath, k, hInv, HInv, ell, workers);

pool = ensureExactParpool(workers);
fprintf('Active parpool workers confirmed: p.NumWorkers=%d (requested %d).\n', ...
    pool.NumWorkers, workers);
checkpoint('parpool', matPath, k, hInv, HInv, ell, workers);

tMesh = tic;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], h);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], H);
meshSeconds = toc(tMesh);
fprintf('Mesh built: N=%d, NT=%d, Nc=%d, NTH=%d, seconds=%.3f\n', ...
    size(node, 1), size(elem, 1), size(nodeH, 1), size(elemH, 1), meshSeconds);
checkpoint('mesh', matPath, k, hInv, HInv, ell, workers, ...
    size(node, 1), size(elem, 1), size(nodeH, 1), size(elemH, 1), meshSeconds);

opts = struct('oversampling', ell, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', true);

tLod = tic;
fprintf('Calling buildLODHelmholtz2D with useParfor=true ...\n');
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, node, elem, bdFlag, ...
    k, 0, 0, opts);
lodSeconds = toc(tLod);

summary = struct();
summary.stage = 'lod_complete';
summary.k = k;
summary.hInv = hInv;
summary.HInv = HInv;
summary.ell = ell;
summary.workers = workers;
summary.N = size(node, 1);
summary.NT = size(elem, 1);
summary.Nc = size(nodeH, 1);
summary.NTH = size(elemH, 1);
summary.lodSeconds = lodSeconds;
summary.nnzTrial = nnz(lod.basis.trial);
summary.nnzTest = nnz(lod.basis.test);
summary.nnzAH = nnz(lod.system.AH);
summary.patchDofMax = max([lod.patch.stats.patchDof]);
summary.freeDofMax = max([lod.patch.stats.freeDof]);
summary.constraintsMax = max([lod.patch.stats.constraints]);
summary.finishedAt = datestr(now, 31);

save(matPath, 'summary', '-v7.3');
fprintf('LOD complete: seconds=%.3f, nnzTrial=%d, nnzTest=%d, nnzAH=%d\n', ...
    summary.lodSeconds, summary.nnzTrial, summary.nnzTest, summary.nnzAH);
fprintf('Patch max: patchDof=%d, freeDof=%d, constraints=%d\n', ...
    summary.patchDofMax, summary.freeDofMax, summary.constraintsMax);
fprintf('MAT: %s\n', matPath);
fprintf('LOG: %s\n', logPath);
fprintf('Finished: %s\n', summary.finishedAt);
fprintf('========== CIP LOD system parallel probe complete ==========\n');


function pool = ensureExactParpool(workers)
workers = max(1, min(round(workers), feature('numcores')));
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting local parpool with %d workers ...\n', workers);
    pool = parpool('local', workers);
elseif pool.NumWorkers ~= workers
    fprintf('Existing parpool has %d workers; restarting with %d workers ...\n', ...
        pool.NumWorkers, workers);
    delete(pool);
    pool = parpool('local', workers);
else
    fprintf('Using existing parpool with requested %d workers.\n', ...
        pool.NumWorkers);
end
if pool.NumWorkers ~= workers
    error('debug_cip_lod_system_parallel_probe:workers', ...
        'Requested %d workers, but active p.NumWorkers=%d.', ...
        workers, pool.NumWorkers);
end
end


function checkpoint(stage, matPath, k, hInv, HInv, ell, workers, varargin)
state = struct();
state.stage = stage;
state.k = k;
state.hInv = hInv;
state.HInv = HInv;
state.ell = ell;
state.workers = workers;
state.time = datestr(now, 31);
if numel(varargin) >= 5
    state.N = varargin{1};
    state.NT = varargin{2};
    state.Nc = varargin{3};
    state.NTH = varargin{4};
    state.meshSeconds = varargin{5};
end
save(matPath, 'state', '-v7.3');
fprintf('Checkpoint: %s at %s\n', stage, state.time);
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
