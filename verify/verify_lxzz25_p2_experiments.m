% VERIFY_LXZZ25_P2_EXPERIMENTS  P2 fine-space LXZZ/LOD-DDM non-paper sweep.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

cfg = p2Config();
outDir = fullfile(repoRoot, 'verify', 'lxzz25_hybrid_lod');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cases = p2Cases(cfg);
results = repmat(emptyResult(), numel(cases), 1);

fprintf('========== LXZZ25 P2 fine-space LOD-DDM experiments ==========\n\n');
fprintf('Memory gate %.1f GB, permission gate %.1f GB, time limit %.0fs.\n', ...
    cfg.memoryLimitGB, cfg.permissionMemoryGB, cfg.timeLimitS);

for i = 1:numel(cases)
    c = cases(i);
    est = estimateCase(c, cfg);
    fprintf('[%d/%d] k=%g %s h=1/%d H=1/%d P2 dofs=%d est=%.2f GB\n', ...
        i, numel(cases), c.k, c.variant, c.hInv, c.HInv, est.p2Dof, est.totalGB);

    r = resultFromCase(c, est);
    if est.totalGB > cfg.memoryLimitGB
        r.status = 'blocked_memory_gt_300GB';
        r.notes = sprintf('Estimate %.2f GB exceeds %.2f GB memory limit.', ...
            est.totalGB, cfg.memoryLimitGB);
        results(i) = r;
        writeOutputs(results, outDir, cfg);
        continue;
    end
    if est.totalGB > cfg.permissionMemoryGB && ~cfg.allowPermissionRows
        r.status = 'needs_permission_gt_200GB';
        r.notes = sprintf('Estimate %.2f GB exceeds %.2f GB permission gate.', ...
            est.totalGB, cfg.permissionMemoryGB);
        results(i) = r;
        writeOutputs(results, outDir, cfg);
        continue;
    end
    if (~cfg.runAllPermitted) && est.p2Dof > cfg.maxRunDof
        r.status = 'queued_runtime_cap';
        r.notes = sprintf('Below memory gate but above interactive cap: P2 dofs=%d.', est.p2Dof);
        results(i) = r;
        writeOutputs(results, outDir, cfg);
        continue;
    end
    if ~cfg.runEnabled
        r.status = 'planned';
        r.notes = 'Run disabled by LXZZ25_P2_RUN=0.';
        results(i) = r;
        writeOutputs(results, outDir, cfg);
        continue;
    end

    try
        runStart = tic;
        out = runOneCase(c, cfg);
        if toc(runStart) > cfg.timeLimitS
            error('verify_lxzz25_p2:timeLimit', ...
                'Case exceeded %.0f seconds after completion.', cfg.timeLimitS);
        end
        r.gmresIter = out.gmresIter;
        r.flag = out.flag;
        r.finalRelres = out.finalRelres;
        r.resvecLength = out.resvecLength;
        r.coarseSetupS = out.coarseSetupS;
        r.localSetupS = out.localSetupS;
        r.solveS = out.solveS;
        r.localMode = out.localMode;
        r.localStoredLuEstimateGB = out.localStoredLuEstimateGB;
        r.localDofMin = out.localDofMin;
        r.localDofMax = out.localDofMax;
        r.localDofMean = out.localDofMean;
        r.status = 'ran';
        if out.flag ~= 0
            r.notes = 'GMRES did not reach tolerance before max iteration.';
        end
    catch ME
        r.status = 'error';
        r.notes = ME.message;
    end
    results(i) = r;
    writeOutputs(results, outDir, cfg);
end

fprintf('\nCSV: %s\n', fullfile(outDir, 'lxzz25_p2_results.csv'));
fprintf('Markdown: %s\n', fullfile(outDir, 'lxzz25_p2_results.md'));
fprintf('========== LXZZ25 P2 experiments complete ==========\n');


function cfg = p2Config()
cfg.kValues = envVector('LXZZ25_P2_KVALUES', [16, 32, 64, 128]);
cfg.variants = cellstr(split(string(envString('LXZZ25_P2_VARIANTS', 'dirichlet,impedance')), ',')).';
cfg.memoryLimitGB = envNumber('LXZZ25_P2_MEMORY_LIMIT_GB', 300);
cfg.permissionMemoryGB = envNumber('LXZZ25_P2_PERMISSION_MEMORY_GB', 200);
cfg.allowPermissionRows = logical(envNumber('LXZZ25_P2_ALLOW_GT_200', 0));
cfg.maxRunDof = envNumber('LXZZ25_P2_MAX_RUN_DOF', 100000);
cfg.runAllPermitted = logical(envNumber('LXZZ25_P2_RUN_ALL_PERMITTED', 0));
cfg.runEnabled = logical(envNumber('LXZZ25_P2_RUN', 1));
cfg.useParfor = logical(envNumber('LXZZ25_P2_PARFOR', 0));
cfg.tol = envNumber('LXZZ25_P2_TOL', 1e-6);
cfg.maxit = envNumber('LXZZ25_P2_MAXIT', 100);
cfg.timeLimitS = envNumber('LXZZ25_P2_TIME_LIMIT_S', 7200);
cfg.m = envNumber('LXZZ25_P2_M', 2);
cfg.lodSolverMode = envString('LXZZ25_P2_LOD_SOLVER_MODE', 'direct');
cfg.solverMode = envString('LXZZ25_P2_SOLVER_MODE', 'adaptive');
cfg.localStoredLuLimitGB = envNumber('LXZZ25_P2_LOCAL_LU_LIMIT_GB', 200);
cfg.luFillConstant = envNumber('LXZZ25_P2_LU_FILL', 40);
cfg.gmresBasisLength = envNumber('LXZZ25_P2_GMRES_BASIS', 103);
end


function cases = p2Cases(cfg)
cases = repmat(emptyCase(), 0, 1);
for ik = 1:numel(cfg.kValues)
    k = cfg.kValues(ik);
    hInv = alignFineInv(ceil(k^(3/2)), k);
    for iv = 1:numel(cfg.variants)
        c = emptyCase();
        c.k = k;
        c.variant = strtrim(cfg.variants{iv});
        c.degree = 2;
        c.hInv = hInv;
        c.HInv = k;
        c.m = cfg.m;
        cases(end+1) = c; %#ok<AGROW>
    end
end
end


function c = emptyCase()
c = struct('k', NaN, 'variant', '', 'degree', 2, 'hInv', NaN, 'HInv', NaN, 'm', 2);
end


function r = emptyResult()
r = struct('k', NaN, 'variant', '', 'degree', NaN, 'hInv', NaN, 'HInv', NaN, ...
    'm', NaN, 'p2Dof', NaN, 'nSub', NaN, 'estimateGB', NaN, ...
    'localMode', '', 'localStoredLuEstimateGB', NaN, 'localDofMin', NaN, ...
    'localDofMax', NaN, 'localDofMean', NaN, 'gmresIter', NaN, 'flag', NaN, ...
    'finalRelres', NaN, 'resvecLength', NaN, 'coarseSetupS', NaN, ...
    'localSetupS', NaN, 'solveS', NaN, 'status', '', 'notes', '');
end


function est = estimateCase(c, cfg)
n = c.hInv;
p2Dof = (2 * n + 1)^2;
nt = 2 * n^2;
nSub = (c.HInv + 1)^2;
supportSteps = max(1, ceil(2 * n / c.HInv));
localDof = (2 * supportSteps + 1)^2;
globalBytes = 3 * 25 * p2Dof * 16;
topologyBytes = 48 * nt;
gmresBytes = 16 * p2Dof * (cfg.gmresBasisLength + 3);
localMatrixBytes = nSub * 3 * 25 * localDof * 16;
localLuBytes = nSub * 16 * cfg.luFillConstant * localDof * log2(max(localDof, 2));
coarseBytes = 3 * 25 * (c.HInv + 1)^2 * 16;
lodBasisBytes = 2 * nnzEstimateEmbeddedLod(p2Dof, (c.HInv + 1)^2, c.m) * 16;
totalBytes = globalBytes + topologyBytes + gmresBytes + localMatrixBytes + ...
    localLuBytes + coarseBytes + lodBasisBytes;

est = struct();
est.p2Dof = p2Dof;
est.nSub = nSub;
est.localDof = localDof;
est.localLuGB = localLuBytes / 1024^3;
est.totalGB = totalBytes / 1024^3;
end


function nz = nnzEstimateEmbeddedLod(nFine, nCoarse, m)
patchFactor = min(nFine, max(10, (2 * m + 3)^2) * nCoarse);
nz = patchFactor;
end


function out = runOneCase(c, cfg)
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1 / c.hInv);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1 / c.HInv);
parts = coarseHatPartition2D(node, elem, bdFlag, 1 / c.HInv);

lodOpts = struct('oversampling', c.m, 'solveCoarse', false, ...
    'solverMode', cfg.lodSolverMode, 'useParfor', cfg.useParfor);
preOpts = struct('degree', 2, 'variant', c.variant, 'coarseType', 'lod', ...
    'lodOptions', lodOpts, 'solverMode', cfg.solverMode, ...
    'localStoredLuLimitGB', cfg.localStoredLuLimitGB, ...
    'localLuFillConstant', cfg.luFillConstant, 'useParfor', cfg.useParfor, ...
    'adjointType', 'reference');

pre = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, c.k, ...
    parts, nodeH, elemH, bdH, preOpts);
b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.k, 2);
tSolve = tic;
[~, flag, relres, iter, resvec] = gmres(pre.A, b, [], cfg.tol, cfg.maxit, @pre.applyResidual);
solveS = toc(tSolve);

out = struct();
out.flag = flag;
out.finalRelres = relres;
out.gmresIter = gmresIterationCount(iter);
out.resvecLength = numel(resvec);
out.coarseSetupS = pre.timing.coarseSetup;
out.localSetupS = pre.timing.localSetup;
out.solveS = solveS;
out.localMode = pre.local.solverModeEffective;
out.localStoredLuEstimateGB = pre.local.estimatedStoredLuGB;
out.localDofMin = pre.local.localDofMin;
out.localDofMax = pre.local.localDofMax;
out.localDofMean = pre.local.localDofMean;
end


function r = resultFromCase(c, est)
r = emptyResult();
r.k = c.k;
r.variant = c.variant;
r.degree = c.degree;
r.hInv = c.hInv;
r.HInv = c.HInv;
r.m = c.m;
r.p2Dof = est.p2Dof;
r.nSub = est.nSub;
r.estimateGB = est.totalGB;
r.localStoredLuEstimateGB = est.localLuGB;
r.localMode = 'estimated';
r.status = 'pending';
end


function writeOutputs(results, outDir, cfg)
csvPath = fullfile(outDir, 'lxzz25_p2_results.csv');
mdPath = fullfile(outDir, 'lxzz25_p2_results.md');

fid = fopen(csvPath, 'w');
cleanup = onCleanup(@() fclose(fid));
names = fieldnames(results);
fprintf(fid, '%s\n', strjoin(names, ','));
for i = 1:numel(results)
    vals = cell(1, numel(names));
    for j = 1:numel(names)
        v = results(i).(names{j});
        if isnumeric(v)
            vals{j} = num2str(v, 16);
        else
            vals{j} = csvEscape(v);
        end
    end
    fprintf(fid, '%s\n', strjoin(vals, ','));
end
clear cleanup;

fid = fopen(mdPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: LXZZ25-inspired P2 fine-space LOD-DDM experiments\n');
fprintf(fid, 'Created: 2026-06-03\n');
fprintf(fid, 'Updated: 2026-06-03\n');
fprintf(fid, 'Verification entry point: `verify/verify_lxzz25_p2_experiments.m`\n');
fprintf(fid, 'Main utilities: `twoLevelHybridSchwarzHelmholtzLOD2D`, `assemblePlaneWaveBoundaryLoad2D`, `coarseHatPartition2D`, `prolongate_P1_P2`, MATLAB `gmres`\n\n');
fprintf(fid, 'Configuration: P2 fine space, P1 LOD basis embedded by $E_{21}$, $m=%g$, Euclidean reference adjoint, memory gate %.1f GB, permission gate %.1f GB, time limit %.0f seconds.\n\n', ...
    cfg.m, cfg.memoryLimitGB, cfg.permissionMemoryGB, cfg.timeLimitS);
fprintf(fid, '| k | variant | degree | h^{-1} | H^{-1} | m | P2 dofs | estimate GB | local mode | LU estimate GB | GMRES it | relres | status | notes |\n');
fprintf(fid, '|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %g | %s | %g | %g | %g | %g | %g | %.2f | %s | %.2f | %s | %s | %s | %s |\n', ...
        r.k, r.variant, r.degree, r.hInv, r.HInv, r.m, r.p2Dof, ...
        r.estimateGB, r.localMode, r.localStoredLuEstimateGB, ...
        mdNumber(r.gmresIter), mdNumber(r.finalRelres), r.status, mdEscape(r.notes));
end
clear cleanup;
end


function hInv = alignFineInv(raw, divisors)
hInv = max(1, ceil(raw));
for d = divisors(:).'
    if d <= 0, continue; end
    hInv = ceil(hInv / d) * d;
end
end


function n = gmresIterationCount(iter)
if numel(iter) == 2
    n = iter(2);
else
    n = iter;
end
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


function s = csvEscape(v)
s = char(string(v));
if contains(s, ',') || contains(s, '"')
    s = ['"', strrep(s, '"', '""'), '"'];
end
end


function s = mdNumber(v)
if isnumeric(v) && isfinite(v)
    s = sprintf('%.4g', v);
else
    s = '-';
end
end


function s = mdEscape(v)
s = char(string(v));
s = strrep(s, '|', '/');
end
