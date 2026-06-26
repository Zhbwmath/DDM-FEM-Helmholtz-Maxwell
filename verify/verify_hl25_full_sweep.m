% VERIFY_HL25_FULL_SWEEP  Checkpointed Hu-Li/LXZZ sweep for kappa/pi 16--128.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

outDir = fullfile(repoRoot, 'tasks', 'HL25_Helmholtz_harmonic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

oldEnv = captureEnv();
cleanup = onCleanup(@() restoreEnv(oldEnv));

fullKappaPi = envString('HL25_FULL_KAPPA_PI', '16,32,64,128');
fullBeta = envNumberList('HL25_FULL_BETA', [0.5, 0.6, 0.7]);

setDefaultEnv('HL25_CROSS_KAPPA_PI', fullKappaPi);
setDefaultEnv('HL25_CROSS_OUTPUT_STEM', 'full_sweep_lxzz_cross_results');
setDefaultEnv('HL25_CROSS_FIGURE_STEM', 'full_sweep_lxzz_cross_iterations');
setDefaultEnv('HL25_CROSS_RUN', '1');
setDefaultEnv('HL25_CROSS_RESUME', '1');
setDefaultEnv('HL25_CROSS_RUN_ALL_PERMITTED', '1');
setDefaultEnv('HL25_CROSS_MAX_RUN_DOF', 'Inf');
setDefaultEnv('HL25_CROSS_PARFOR', '1');
setDefaultEnv('HL25_CROSS_PARFOR_WORKERS', ...
    envString('HL25_FULL_PARFOR_WORKERS', '12'));
setDefaultEnv('HL25_CROSS_SOLVER_MODE', 'adaptive');
setDefaultEnv('HL25_CROSS_COARSE_SOLVER_MODE', 'adaptive');
setDefaultEnv('HL25_CROSS_RANK_METHOD', 'none');
setDefaultEnv('HL25_CROSS_PERMISSION_MEMORY_GB', '200');
setDefaultEnv('HL25_CROSS_HARD_MEMORY_GB', '500');
setDefaultEnv('HL25_CROSS_TIME_LIMIT_S', '86400');

fprintf('========== Hu-Li / LXZZ25 full sweep ==========\n');
fprintf('kappa/pi: %s\n', getenv('HL25_CROSS_KAPPA_PI'));
fprintf('beta values: %s\n', num2str(fullBeta));
fprintf('output stem: %s\n\n', getenv('HL25_CROSS_OUTPUT_STEM'));

for i = 1:numel(fullBeta)
    setenv('HL25_CROSS_BETA', sprintf('%.16g', fullBeta(i)));
    fprintf('\n---------- beta = %.16g ----------\n', fullBeta(i));
    run(fullfile(repoRoot, 'verify', 'verify_hl25_lxzz_cross_study.m'));
end

fprintf('\nCSV: %s\n', fullfile(outDir, ...
    [getenv('HL25_CROSS_OUTPUT_STEM'), '.csv']));
fprintf('Markdown: %s\n', fullfile(outDir, ...
    [getenv('HL25_CROSS_OUTPUT_STEM'), '.md']));
fprintf('========== Hu-Li / LXZZ25 full sweep complete ==========\n');


function oldEnv = captureEnv()
names = {'HL25_CROSS_KAPPA_PI', 'HL25_CROSS_OUTPUT_STEM', ...
    'HL25_CROSS_FIGURE_STEM', 'HL25_CROSS_RUN', ...
    'HL25_CROSS_RESUME', 'HL25_CROSS_RUN_ALL_PERMITTED', ...
    'HL25_CROSS_MAX_RUN_DOF', 'HL25_CROSS_PARFOR', ...
    'HL25_CROSS_PARFOR_WORKERS', ...
    'HL25_CROSS_SOLVER_MODE', 'HL25_CROSS_COARSE_SOLVER_MODE', ...
    'HL25_CROSS_RANK_METHOD', 'HL25_CROSS_PERMISSION_MEMORY_GB', ...
    'HL25_CROSS_HARD_MEMORY_GB', 'HL25_CROSS_TIME_LIMIT_S', ...
    'HL25_CROSS_BETA'};
oldEnv.names = names;
oldEnv.values = cell(size(names));
oldEnv.wasSet = false(size(names));
for i = 1:numel(names)
    value = getenv(names{i});
    oldEnv.values{i} = value;
    oldEnv.wasSet(i) = ~isempty(value);
end
end


function restoreEnv(oldEnv)
for i = 1:numel(oldEnv.names)
    if oldEnv.wasSet(i)
        setenv(oldEnv.names{i}, oldEnv.values{i});
    else
        setenv(oldEnv.names{i}, '');
    end
end
end


function setDefaultEnv(name, value)
if isempty(getenv(name))
    setenv(name, value);
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
