% DEBUG_CIP_LXZZ_LOCAL_APPLY_MODES  Compare full and compact local apply modes.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

cfg = localConfig();
outDir = fullfile(repoRoot, 'tasks', 'Frmwrk_hybrd2lvl_DDM');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('========== CIP/LXZZ local apply mode probe ==========\n');
fprintf('k=%g degree=P%d h=1/%d localH=1/%d workers=%d rhs=%d\n', ...
    cfg.k, cfg.degree, cfg.hInv, cfg.localHInv, cfg.workers, cfg.nRhs);

ensurePool(cfg.workers);

[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1 / cfg.hInv);
parts = coarseHatPartition2D(node, elem, bdFlag, 1 / cfg.localHInv);
fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, cfg.k, ...
    struct('degree', cfg.degree, 'cacheEnergySolver', false));

rng(11);
r = randn(fine.N, cfg.nRhs) + 1i * randn(fine.N, cfg.nRhs);

records = repmat(emptyRecord(), 3, 1);
records(1) = runMode('serial', fine, parts, r, cfg, false, 'full');
records(2) = runMode('parfor-full', fine, parts, r, cfg, true, 'full');
records(3) = runMode('parfor-compact', fine, parts, r, cfg, true, 'compact');

relCompact = norm(records(2).y - records(3).y, 'fro') / ...
    max(1, norm(records(2).y, 'fro'));
relSerial = norm(records(1).y - records(2).y, 'fro') / ...
    max(1, norm(records(1).y, 'fro'));
fprintf('serial/full relative difference: %.3e\n', relSerial);
fprintf('full/compact relative difference: %.3e\n', relCompact);

writeReport(outDir, records, cfg, relSerial, relCompact, numel(parts), fine.N);
fprintf('Report: %s\n', fullfile(outDir, 'local_apply_mode_probe.md'));
fprintf('========== local apply mode probe complete ==========\n');


function cfg = localConfig()
cfg.k = envNumber('CIP_LXZZ_APPLY_PROBE_K', 8);
cfg.degree = envNumber('CIP_LXZZ_APPLY_PROBE_DEGREE', 1);
cfg.hInv = envNumber('CIP_LXZZ_APPLY_PROBE_HINV', 16);
cfg.localHInv = envNumber('CIP_LXZZ_APPLY_PROBE_LOCAL_HINV', 4);
cfg.workers = envNumber('CIP_LXZZ_APPLY_PROBE_WORKERS', 2);
cfg.nRhs = envNumber('CIP_LXZZ_APPLY_PROBE_NRHS', 2);
cfg.solverMode = envString('CIP_LXZZ_APPLY_PROBE_SOLVER_MODE', 'direct');
cfg.localStorage = envString('CIP_LXZZ_APPLY_PROBE_LOCAL_STORAGE', 'matrix');
cfg.repeat = envNumber('CIP_LXZZ_APPLY_PROBE_REPEAT', 3);
end


function record = runMode(label, fine, parts, r, cfg, applyParfor, applyMode)
opts = struct('solverMode', cfg.solverMode, ...
    'localStorage', cfg.localStorage, ...
    'useParfor', applyParfor, ...
    'setupParfor', false, ...
    'applyParfor', applyParfor, ...
    'applyMode', applyMode);
local = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, ...
    'dirichlet', opts);
local.applyInverse(r);
memBefore = matlabMemoryGB();
times = zeros(cfg.repeat, 1);
for j = 1:cfg.repeat
    t = tic;
    y = local.applyInverse(r);
    times(j) = toc(t);
end
memAfter = matlabMemoryGB();
record = emptyRecord();
record.label = label;
record.applyMode = applyMode;
record.applyParfor = applyParfor;
record.secondsMedian = median(times);
record.secondsMin = min(times);
record.memBeforeGB = memBefore;
record.memAfterGB = memAfter;
record.y = y;
record.localDofMax = local.info.localDofMax;
record.localDofMean = local.info.localDofMean;
end


function record = emptyRecord()
record = struct('label', '', 'applyMode', '', 'applyParfor', false, ...
    'secondsMedian', NaN, 'secondsMin', NaN, 'memBeforeGB', NaN, ...
    'memAfterGB', NaN, 'localDofMax', NaN, 'localDofMean', NaN, ...
    'y', []);
end


function writeReport(outDir, records, cfg, relSerial, relCompact, nSub, nGlobal)
path = fullfile(outDir, 'local_apply_mode_probe.md');
pool = gcp('nocreate');
if isempty(pool)
    nWorkers = 1;
else
    nWorkers = pool.NumWorkers;
end
nBlocks = min(nSub, max(1, 2 * nWorkers));
fullOutputGB = 16 * nGlobal * cfg.nRhs * nBlocks / 2^30;

fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Local apply mode timing probe for CIP/LXZZ local solvers.\n');
fprintf(fid, 'Created: 2026-07-03\n');
fprintf(fid, 'Updated: 2026-07-03\n');
fprintf(fid, 'Verification entry point: `debug/debug_cip_lxzz_local_apply_modes.m`\n');
fprintf(fid, 'Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, MATLAB `parfor`\n\n');
fprintf(fid, '# Local Apply Mode Probe\n\n');
fprintf(fid, 'Settings: $k=%g$, P%d, $h^{-1}=%d$, local $H^{-1}=%d$, workers `%d`, right-hand sides `%d`, repeats `%d`.\n\n', ...
    cfg.k, cfg.degree, cfg.hInv, cfg.localHInv, nWorkers, cfg.nRhs, cfg.repeat);
fprintf(fid, 'Estimated full-vector worker-output memory: %.4f GB for `%d` blocks and `%d` global DOFs.\n\n', ...
    fullOutputGB, nBlocks, nGlobal);
fprintf(fid, '| mode | parfor | apply mode | median s | min s | MATLAB memory before GB | after GB | local dof max | local dof mean |\n');
fprintf(fid, '|---|---:|---|---:|---:|---:|---:|---:|---:|\n');
for i = 1:numel(records)
    r = records(i);
    fprintf(fid, '| %s | %d | %s | %.4g | %.4g | %.4g | %.4g | %.4g | %.4g |\n', ...
        r.label, r.applyParfor, r.applyMode, r.secondsMedian, ...
        r.secondsMin, r.memBeforeGB, r.memAfterGB, ...
        r.localDofMax, r.localDofMean);
end
fprintf(fid, '\nRelative differences: serial/full %.3e; full/compact %.3e.\n', ...
    relSerial, relCompact);
fprintf(fid, '\nInterpretation rule: keep `full` for cases where the worker-output estimate is below the configured threshold and timing is comparable; use `compact` when the full-vector output estimate crosses the memory limit.\n');
clear cleanup
end


function ensurePool(nWorkers)
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', nWorkers);
elseif pool.NumWorkers ~= nWorkers
    delete(pool);
    parpool('local', nWorkers);
end
end


function gb = matlabMemoryGB()
try
    m = memory;
    gb = m.MemUsedMATLAB / 2^30;
catch
    gb = NaN;
end
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


function s = envString(name, defaultValue)
s = getenv(name);
if isempty(s), s = defaultValue; end
end
