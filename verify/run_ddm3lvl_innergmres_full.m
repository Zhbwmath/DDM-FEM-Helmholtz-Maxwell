% RUN_DDM3LVL_INNERGMRES_FULL  Execute the full source LOD-DDM inner-GMRES sweep.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'DDM3lvl_LOD_Helmholtz');
if ~exist(outDir, 'dir'), mkdir(outDir); end

diaryPath = fullfile(outDir, 'innergmres_lu_full_run.diary.log');
if exist(diaryPath, 'file'), delete(diaryPath); end
diary(diaryPath);
cleanup = onCleanup(@() diary('off'));

setenv('DDM3LVL_RUN', '1');
setenv('DDM3LVL_K_VALUES', '16 32 64');
setenv('DDM3LVL_H0_INV_VALUES', '1 2 4 8');
setenv('DDM3LVL_M_VALUES', '1 2 3 4');
setenv('DDM3LVL_MAX_RUN_DOF', '1000000');
setenv('DDM3LVL_MAXIT', '80');
setenv('DDM3LVL_OUTPUT_SUFFIX', 'innergmres_lu_full');
setenv('DDM3LVL_COARSE_SOLVE_MODE', 'innerGmres');
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
        error('run_ddm3lvl_innergmres_full:workers', ...
            'DDM3LVL_WORKERS must be a positive integer.');
    end
    workers = floor(workers);
end
setenv('DDM3LVL_PARPOOL_WORKERS', num2str(workers));
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', workers);
end
if pool.NumWorkers ~= workers
    error('run_ddm3lvl_innergmres_full:poolSize', ...
        'Expected %d workers, found %d.', workers, pool.NumWorkers);
end

fprintf('Started full inner-GMRES source LOD-DDM sweep at %s.\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
run(fullfile(repoRoot, 'verify', 'verify_ddm3lvl_lod_helmholtz_experiments.m'));
fprintf('Finished full inner-GMRES source LOD-DDM sweep at %s.\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
