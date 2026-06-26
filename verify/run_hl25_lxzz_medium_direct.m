% RUN_HL25_LXZZ_MEDIUM_DIRECT  Direct MATLAB entry point for the medium Hu-Li/LXZZ sweep.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
cd(repoRoot);
addpath(genpath(repoRoot));

setDefaultEnv('HL25_LXZZ_MEDIUM_KVALUES', '16,32,64,128');
setDefaultEnv('HL25_LXZZ_MEDIUM_COARSE', 'economic,spectral');
setDefaultEnv('HL25_LXZZ_MEDIUM_VARIANTS', 'dirichlet,impedance');
setDefaultEnv('HL25_LXZZ_MEDIUM_RESUME', '1');
setDefaultEnv('HL25_LXZZ_MEDIUM_RUN', '1');
setDefaultEnv('HL25_LXZZ_MEDIUM_MAX_RUN_DOF', '3000000');
setDefaultEnv('HL25_LXZZ_MEDIUM_MAX_COARSE', '300000');
setDefaultEnv('HL25_LXZZ_MEDIUM_PERMISSION_GB', '200');
setDefaultEnv('HL25_LXZZ_MEDIUM_HARD_GB', '500');
setDefaultEnv('HL25_LXZZ_MEDIUM_OUTPUT_STEM', ...
    'lxzz_hl_coarse_medium_direct_results');
setDefaultEnv('HL25_LXZZ_MEDIUM_PARFOR', '1');
setDefaultEnv('HL25_LXZZ_MEDIUM_WORKERS', '8');
setDefaultEnv('HL25_LXZZ_MEDIUM_PARPOOL_PROFILE', 'zhbw_cluster');

configurePoolFromEnv();

run(fullfile(repoRoot, 'verify', 'verify_hl25_lxzz_hybrid_medium.m'));


function setDefaultEnv(name, value)
if isempty(getenv(name))
    setenv(name, value);
end
end


function configurePoolFromEnv()
useParfor = str2double(getenv('HL25_LXZZ_MEDIUM_PARFOR')) ~= 0;
if ~useParfor
    return;
end
workerCount = str2double(getenv('HL25_LXZZ_MEDIUM_WORKERS'));
if ~isfinite(workerCount) || workerCount <= 0
    return;
end
workerCount = max(1, round(workerCount));
profile = getenv('HL25_LXZZ_MEDIUM_PARPOOL_PROFILE');
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers == workerCount
    return;
end
if ~isempty(pool)
    delete(pool);
end
if isempty(profile)
    parpool(workerCount);
else
    parpool(profile, workerCount);
end
end
