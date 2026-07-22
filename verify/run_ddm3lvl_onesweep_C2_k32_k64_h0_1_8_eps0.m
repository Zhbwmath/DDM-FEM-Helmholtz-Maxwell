% RUN_DDM3LVL_ONESWEEP_C2_K32_K64_H0_1_8_EPS0  Rerun requested C=2 sample rows.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'DDM3lvl_LOD_Helmholtz');
if ~exist(outDir, 'dir'), mkdir(outDir); end

diaryPath = fullfile(outDir, 'onesweep_C2_k32_k64_h0_1_8_eps0_run.diary.log');
if exist(diaryPath, 'file'), delete(diaryPath); end
diary(diaryPath);
cleanup = onCleanup(@() diary('off'));

setenv('DDM3LVL_RUN', '1');
setenv('DDM3LVL_K_VALUES', '32 64');
setenv('DDM3LVL_H0_INV_VALUES', '8');
setenv('DDM3LVL_M_VALUES', '4');
setenv('DDM3LVL_INCLUDE_LOG_M', '0');
setenv('DDM3LVL_EPSILON_LABELS', 'zero');
setenv('DDM3LVL_MAX_RUN_DOF', '1000000');
setenv('DDM3LVL_MAXIT', '80');
setenv('DDM3LVL_OUTPUT_SUFFIX', 'onesweep_C2_k32_k64_h0_1_8_eps0');
setenv('DDM3LVL_TILDE_BUFFER_LAYERS', '2');
setenv('DDM3LVL_COARSE_SOLVE_MODE', 'oneSweep');
setenv('DDM3LVL_COARSE_FACTOR_MODE', 'lu');
setenv('DDM3LVL_COARSE_INNER_TOL', '1e-6');
setenv('DDM3LVL_COARSE_INNER_MAXIT', '80');
setenv('DDM3LVL_PARALLEL_LOCAL_SETUP', '1');
setenv('DDM3LVL_COMPARE_LOCAL_BASIS', '0');

workerText = getenv('DDM3LVL_WORKERS');
if isempty(workerText)
    workers = min(max(1, feature('numcores')), 12);
else
    workers = str2double(workerText);
    if ~isfinite(workers) || workers < 1
        error('run_ddm3lvl_onesweep_C2_k32_k64_h0_1_8_eps0:workers', ...
            'DDM3LVL_WORKERS must be a positive integer.');
    end
    workers = floor(workers);
end
setenv('DDM3LVL_PARPOOL_WORKERS', num2str(workers));
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', workers);
elseif pool.NumWorkers ~= workers
    delete(pool);
    pool = parpool('local', workers);
end
if pool.NumWorkers ~= workers
    error('run_ddm3lvl_onesweep_C2_k32_k64_h0_1_8_eps0:poolSize', ...
        'Expected %d workers, found %d.', workers, pool.NumWorkers);
end

fprintf('Started C=2 one-shot source LOD-DDM H0=1/8 eps=0 sample at %s.\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
run(fullfile(repoRoot, 'verify', 'verify_ddm3lvl_lod_helmholtz_experiments.m'));
fprintf('Finished C=2 one-shot source LOD-DDM H0=1/8 eps=0 sample at %s.\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
