% VERIFY_CIP_LXZZ_LOD_MEDIUM  CIP fine form with LXZZ LOD coarse space.

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

fprintf('========== CIP fine form with LXZZ-type LOD coarse space ==========\n\n');
fprintf('k values: %s\n', mat2str(cfg.kValues));
fprintf('degree P%d, C_h=%g, max run dof=%g, run=%d, force-run=%d, parfor=%d\n\n', ...
    cfg.degree, cfg.Ch, cfg.maxRunDof, cfg.runEnabled, cfg.forceRun, ...
    cfg.useParfor);

ensureParpool(cfg);

if cfg.runSmoke
    runSmallIdentitySmoke(cfg);
end

for i = 1:numel(cases)
    c = finalizeCase(cases(i), cfg);
    if cfg.resume && ~cfg.forceRun && caseComplete(results, c)
        fprintf('[%02d/%02d] k=%g variant=%s: checkpoint complete, skipping.\n', ...
            i, numel(cases), c.k, c.variant);
        continue;
    end

    est = estimateCase(c, cfg);
    row = baseRow(c, est, cfg);
    fprintf('[%02d/%02d] k=%g variant=%s: N=%d, H=1/%d, m=%d, %.2f GB ... ', ...
        i, numel(cases), c.k, c.variant, est.N, c.HInv, c.oversampling, ...
        est.totalGB);

    if cfg.forceRun
        row.status = 'running';
        row.notes = appendNote(row.notes, sprintf( ...
            'force-run started; estimate %.2f GB.', est.totalGB));
        results = replaceRow(results, row);
        writeOutputs(outDir, results, cfg);
        try
            row = runCase(row, c, cfg);
            row.notes = appendNote(row.notes, sprintf( ...
                'force-run bypassed gates; estimate %.2f GB.', est.totalGB));
            fprintf('FORCE-RAN: %d iterations, relres %.3e\n', ...
                row.gmresIterations, row.finalRelres);
        catch ME
            row.status = 'failed';
            row.notes = errorNote(ME);
            fprintf('FAILED: %s\n%s\n', ME.message, errorReport(ME));
        end
    elseif est.totalGB > cfg.hardMemoryGB
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
        row.notes = sprintf('Above medium run cap N=%d.', cfg.maxRunDof);
        fprintf('QUEUED DOF\n');
    elseif ~cfg.runEnabled
        row.status = 'estimated_only';
        row.notes = 'Execution disabled by CIP_LXZZ_LOD_RUN=0.';
        fprintf('ESTIMATED\n');
    else
        try
            row = runCase(row, c, cfg);
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
fprintf('========== CIP/LXZZ LOD medium verification complete ==========\n');


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
cfg.kValues = envVector('CIP_LXZZ_LOD_KVALUES', [16, 32, 64, 128]);
cfg.variants = envList('CIP_LXZZ_LOD_VARIANTS', {'dirichlet', 'impedance'});
cfg.degree = envNumber('CIP_LXZZ_LOD_DEGREE', 1);
cfg.Ch = envNumber('CIP_LXZZ_LOD_CH', 1);
cfg.gamma = envNumberOrEmpty('CIP_LXZZ_LOD_GAMMA', []);
cfg.cipOptions = struct();
cfg.tol = envNumber('CIP_LXZZ_LOD_TOL', 1e-6);
cfg.maxit = envNumber('CIP_LXZZ_LOD_MAXIT', 100);
cfg.restart = [];
cfg.runEnabled = logical(envNumber('CIP_LXZZ_LOD_RUN', 1));
cfg.resume = logical(envNumber('CIP_LXZZ_LOD_RESUME', 1));
cfg.forceRun = logical(envNumber('CIP_LXZZ_LOD_FORCE_RUN', 0));
cfg.runSmoke = logical(envNumber('CIP_LXZZ_LOD_SMOKE', 1));
cfg.maxRunDof = envNumber('CIP_LXZZ_LOD_MAX_RUN_DOF', 50000);
cfg.permissionMemoryGB = envNumber('CIP_LXZZ_LOD_PERMISSION_GB', 200);
cfg.hardMemoryGB = envNumber('CIP_LXZZ_LOD_HARD_GB', 500);
cfg.allowPermissionRows = logical(envNumber('CIP_LXZZ_LOD_ALLOW_GT_200', 0));
cfg.solverMode = envString('CIP_LXZZ_LOD_SOLVER_MODE', 'lu');
cfg.lodSolverMode = envString('CIP_LXZZ_LOD_LOD_SOLVER', 'direct');
cfg.adjointType = envString('CIP_LXZZ_LOD_ADJOINT', 'energy');
cfg.useParfor = logical(envNumber('CIP_LXZZ_LOD_PARFOR', 1));
cfg.parpoolWorkers = envNumber('CIP_LXZZ_LOD_PARPOOL_WORKERS', feature('numcores'));
cfg.localSetupParfor = logical(envNumber('CIP_LXZZ_LOD_LOCAL_SETUP_PARFOR', 0));
cfg.localApplyParfor = logical(envNumber('CIP_LXZZ_LOD_LOCAL_APPLY_PARFOR', cfg.useParfor));
cfg.localApplyMode = envString('CIP_LXZZ_LOD_LOCAL_APPLY_MODE', 'auto');
cfg.fullVectorApplyLimitGB = envNumber('CIP_LXZZ_LOD_FULL_VECTOR_APPLY_LIMIT_GB', 2);
cfg.lodParfor = logical(envNumber('CIP_LXZZ_LOD_LOD_PARFOR', 0));
cfg.localStorage = envString('CIP_LXZZ_LOD_LOCAL_STORAGE', 'matrix');
cfg.oversamplingOverride = envNumberOrEmpty('CIP_LXZZ_LOD_OVERSAMPLING', []);
cfg.highKThreshold = envNumber('CIP_LXZZ_LOD_HIGHK_THRESHOLD', 128);
cfg.highKOversampling = envNumber('CIP_LXZZ_LOD_HIGHK_OVERSAMPLING', 3);
cfg.localLuFillConstant = envNumber('CIP_LXZZ_LOD_LU_FILL', 40);
cfg.gmresBasisLength = envNumber('CIP_LXZZ_LOD_GMRES_BASIS', 103);
cfg.outputStem = envString('CIP_LXZZ_LOD_OUTPUT_STEM', ...
    'cip_lxzz_lod_medium_results');
end


function ensureParpool(cfg)
needsPool = cfg.useParfor || cfg.localSetupParfor || ...
    cfg.localApplyParfor || cfg.lodParfor;
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
    error('verify_cip_lxzz_lod_medium:parpoolWorkers', ...
        ['Requested %d parpool workers, but the active pool has %d. ' ...
        'Check CIP_LXZZ_LOD_PARPOOL_WORKERS and the local profile.'], ...
        nWorkers, pool.NumWorkers);
end
fprintf('Active parpool workers confirmed: p.NumWorkers=%d (requested %d).\n', ...
    pool.NumWorkers, nWorkers);
end


function cases = localCases(cfg)
cases = repmat(emptyCase(), 0, 1);
for ik = 1:numel(cfg.kValues)
    for iv = 1:numel(cfg.variants)
        c = emptyCase();
        c.k = cfg.kValues(ik);
        c.variant = cfg.variants{iv};
        cases(end+1) = c; %#ok<AGROW>
    end
end
end


function c = emptyCase()
c = struct('k', NaN, 'variant', '', 'degree', NaN, 'hInv', NaN, ...
    'HInv', NaN, 'localHInv', NaN, 'h', NaN, 'H', NaN, ...
    'localH', NaN, 'oversampling', NaN, 'nSubdomains', NaN);
end


function c = finalizeCase(c, cfg)
c.degree = cfg.degree;
c.HInv = max(1, round(c.k));
rawFine = ceil(cfg.Ch * c.k^((2 * c.degree + 1) / (2 * c.degree)));
c.hInv = alignFineInv(rawFine, c.HInv);
c.h = 1 / c.hInv;
c.H = 1 / c.HInv;
c.localH = selectLocalSpacing(c.k, c.variant);
c.localHInv = max(1, round(1 / c.localH));
c.oversampling = selectOversampling(c.k, cfg);
c.nSubdomains = (c.localHInv + 1)^2;
end


function spacing = selectLocalSpacing(k, variant)
switch lower(variant)
    case {'q1', 'dirichlet'}
        spacing = 1 / k;
    case {'q2', 'impedance'}
        spacing = 2 / k;
    otherwise
        error('verify_cip_lxzz_lod_medium:variant', ...
            'Unknown local solver variant "%s".', variant);
end
end


function ell = selectOversampling(k, cfg)
if ~isempty(cfg.oversamplingOverride)
    ell = cfg.oversamplingOverride;
elseif k >= cfg.highKThreshold
    ell = cfg.highKOversampling;
else
    ell = max(1, round(log2(k) - 1));
end
ell = max(1, round(ell));
end


function est = estimateCase(c, cfg)
N = (c.degree * c.hInv + 1)^2;
NT = 2 * c.hInv^2;
Nc = (c.HInv + 1)^2;
supportSteps = max(1, ceil(2 * c.hInv / c.localHInv));
localDof = (c.degree * supportSteps + 1)^2;
globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.gmresBasisLength + 3);
localMatrixBytes = c.nSubdomains * 336 * localDof;
localLuBytes = c.nSubdomains * 16 * cfg.localLuFillConstant * ...
    localDof * log2(max(localDof, 2));
if any(strcmpi(cfg.localStorage, {'matrix', 'lazy', 'deferred'}))
    localStoredBytes = localMatrixBytes;
else
    localStoredBytes = localMatrixBytes + localLuBytes;
end
lodBasisBytes = 2 * 16 * min(N * Nc, Nc * max(10, ...
    (2 * c.oversampling + 3)^2));
coarseBytes = 112 * Nc + 16 * 20 * Nc * log2(max(Nc, 2));
est = struct('N', N, 'NT', NT, 'Nc', Nc, ...
    'supportSteps', supportSteps, 'localDofEstimate', localDof, ...
    'totalGB', (globalBytes + gmresBytes + localStoredBytes + ...
    lodBasisBytes + coarseBytes) / 2^30);
end


function row = runCase(row, c, cfg)
tAll = tic;
fprintf('\n    mesh setup ... ');
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], c.h);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], c.H);
parts = coarseHatPartition2D(node, elem, bdFlag, c.localH);
fprintf('done');

fprintf('\n    CIP fine space ... ');
fineOpts = struct('degree', c.degree, 'gamma', cfg.gamma, ...
    'cipOptions', cfg.cipOptions, ...
    'cacheEnergySolver', strcmpi(cfg.adjointType, 'energy'));
fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, c.k, fineOpts);
fprintf('done');

fprintf('\n    local solver setup (%s storage) ... ', cfg.localStorage);
localOpts = struct('gamma', cfg.gamma, 'cipOptions', cfg.cipOptions, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'setupParfor', cfg.localSetupParfor, ...
    'applyParfor', cfg.localApplyParfor, ...
    'applyMode', cfg.localApplyMode, ...
    'fullVectorApplyLimitGB', cfg.fullVectorApplyLimitGB, ...
    'localLuFillConstant', cfg.localLuFillConstant, ...
    'localStorage', cfg.localStorage);
localSolver = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, ...
    c.variant, localOpts);
fprintf('done');

fprintf('\n    LOD coarse space (m=%d) ... ', c.oversampling);
lodOpts = struct('oversampling', c.oversampling, 'solveCoarse', false, ...
    'solverMode', cfg.lodSolverMode, 'useParfor', cfg.lodParfor);
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, node, elem, bdFlag, ...
    c.k, 0, 0, lodOpts);
fprintf('done');

fprintf('\n    hybrid assembly ... ');
coarseSpace = struct('nativeTrial', lod.basis.trial, ...
    'nativeTest', lod.basis.test, 'embedding', fine.p1ToFine, ...
    'object', lod, 'description', ...
    'normal-FEM LOD basis with CIP fine Galerkin matrix');
preOpts = struct('fineSpace', fine, 'coarseSpace', coarseSpace, ...
    'localSolver', localSolver, 'variant', c.variant, ...
    'solverMode', cfg.solverMode, 'useParfor', cfg.useParfor, ...
    'adjointType', cfg.adjointType, ...
    'localLuFillConstant', cfg.localLuFillConstant);
pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, c.k, ...
    parts, nodeH, elemH, bdH, preOpts);
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
row.cipNonzeros = nnz(fine.cip.matrix);
row.solveS = toc(tSolve);
row.totalS = toc(tAll);
row.setupS = row.totalS - row.solveS;
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
    error('verify_cip_lxzz_lod_medium:hybridIdentity', ...
        'apply(x) and applyResidual(A*x) differ by %.3e.', err);
end
end


function runSmallIdentitySmoke(cfg)
fprintf('Small identity smoke ... ');
c = finalizeCase(struct('k', 8, 'variant', 'dirichlet', ...
    'degree', NaN, 'hInv', NaN, 'HInv', NaN, 'h', NaN, 'H', NaN, ...
    'oversampling', NaN, 'nSubdomains', NaN), cfg);
c.hInv = 16;
c.h = 1 / c.hInv;
c.HInv = 4;
c.H = 1 / c.HInv;
c.oversampling = 1;
row = baseRow(c, estimateCase(c, cfg), cfg);
row = runCase(row, c, cfg);
assert(row.identityError < 1e-9, 'Smoke identity failed.');
fprintf('PASSED (identity %.2e)\n', row.identityError);
end


function row = emptyRow()
row = struct('k', NaN, 'variant', '', 'degree', NaN, ...
    'fineHInv', NaN, 'localHInv', NaN, 'oversampling', NaN, ...
    'ndof', NaN, 'nSubdomains', NaN, 'ncoarse', NaN, ...
    'supportSteps', NaN, 'localDofEstimate', NaN, ...
    'estimateGB', NaN, 'coarseDimension', NaN, 'cipNonzeros', NaN, ...
    'localMode', '', 'localDofMin', NaN, 'localDofMax', NaN, ...
    'localDofMean', NaN, 'gmresIterations', NaN, 'flag', NaN, ...
    'finalRelres', NaN, 'resvecLength', NaN, 'identityError', NaN, ...
    'setupS', NaN, 'solveS', NaN, 'totalS', NaN, ...
    'status', '', 'notes', '');
end


function row = baseRow(c, est, cfg)
row = emptyRow();
row.k = c.k;
row.variant = c.variant;
row.degree = c.degree;
row.fineHInv = c.hInv;
row.localHInv = c.localHInv;
row.oversampling = c.oversampling;
row.ndof = est.N;
row.nSubdomains = c.nSubdomains;
row.ncoarse = est.Nc;
row.supportSteps = est.supportSteps;
row.localDofEstimate = est.localDofEstimate;
row.estimateGB = est.totalGB;
row.notes = gammaText(cfg.gamma);
row.status = 'pending';
end


function text = appendNote(text, note)
if isempty(text)
    text = note;
else
    text = [text, ' ', note];
end
end


function results = loadResults(outDir, cfg)
path = fullfile(outDir, [cfg.outputStem, '.csv']);
if ~exist(path, 'file')
    results = repmat(emptyRow(), 0, 1);
    return;
end
t = readtable(path, 'TextType', 'string');
results = repmat(emptyRow(), height(t), 1);
for i = 1:height(t)
    r = emptyRow();
    names = fieldnames(r);
    for j = 1:numel(names)
        if ~ismember(names{j}, t.Properties.VariableNames), continue; end
        v = t.(names{j})(i);
        if isnumeric(r.(names{j}))
            r.(names{j}) = double(v);
        elseif ismissing(v)
            r.(names{j}) = '';
        else
            r.(names{j}) = char(v);
        end
    end
    results(i) = r;
end
end


function tf = caseComplete(results, c)
tf = false;
for i = 1:numel(results)
    r = results(i);
    if r.k == c.k && strcmp(r.variant, c.variant) && ...
            ismember(r.status, {'ran', 'ran_not_converged', ...
            'queued_runtime_cap', 'requires_permission_gt_200gb', ...
            'blocked_memory_gt_hard_limit'})
        tf = true;
        return;
    end
end
end


function results = replaceRow(results, row)
for i = 1:numel(results)
    r = results(i);
    if r.k == row.k && strcmp(r.variant, row.variant)
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
names = fieldnames(emptyRow());

fid = fopen(csvPath, 'w');
cleanup = onCleanup(@() fclose(fid));
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
clear cleanup

fid = fopen(mdPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: CIP fine form with normal-FEM LXZZ LOD coarse space.\n');
fprintf(fid, 'Created: 2026-06-26\n');
fprintf(fid, 'Updated: 2026-06-26\n');
fprintf(fid, 'Verification entry point: `verify/verify_cip_lxzz_lod_medium.m`\n');
fprintf(fid, 'Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildLODHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# CIP Fine Form With LXZZ-Type LOD Coarse Space\n\n');
fprintf(fid, 'This first research step uses the LXZZ twice-hybrid residual preconditioner. The fine and local matrices use the CIP sesquilinear form; the coarse basis is the normal-FEM LOD Helmholtz basis. The injected coarse matrix is recomputed against the CIP fine matrix. Fine mesh rule: P%d with $h^{-1}=\\operatorname{align}(\\lceil C_h k^{(2p+1)/(2p)}\\rceil,k)$ and $C_h=%g$. LOD uses $H=1/k$ and $m=\\max(1,\\operatorname{round}(\\log_2 k-1))$ below $k=%g$, then fixed $m=%g$ for high-$k$ rows unless `CIP_LXZZ_LOD_OVERSAMPLING` overrides it. LXZZ local partition spacing is $1/k$ for Dirichlet and $2/k$ for impedance. Local storage is `%s`: local sparse matrices are stored and factored inside each GMRES preconditioner apply. Local setup parfor `%d`, local apply parfor `%d`, local apply mode `%s` with full-vector output limit %.2f GB, LOD parfor `%d`. GMRES tolerance %.1e, max iterations %d, medium run cap $N\\le %g$, force-run flag `%d`.\n\n', ...
    cfg.degree, cfg.Ch, cfg.highKThreshold, cfg.highKOversampling, ...
    cfg.localStorage, cfg.localSetupParfor, cfg.localApplyParfor, ...
    cfg.localApplyMode, cfg.fullVectorApplyLimitGB, cfg.lodParfor, ...
    cfg.tol, cfg.maxit, cfg.maxRunDof, cfg.forceRun);
fprintf(fid, '| k | variant | h^{-1} | H^{-1} | m | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |\n');
fprintf(fid, '|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %g | %s | %g | %g | %g | %g | %s | %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
        r.k, r.variant, r.fineHInv, r.localHInv, r.oversampling, ...
        r.ndof, mdNumber(r.coarseDimension), r.estimateGB, ...
        mdNumber(r.cipNonzeros), r.localMode, mdNumber(r.localDofMax), ...
        mdNumber(r.identityError), mdNumber(r.gmresIterations), ...
        mdNumber(r.finalRelres), mdNumber(r.setupS), mdNumber(r.solveS), ...
        r.status, mdEscape(r.notes));
end
clear cleanup
end


function hInv = alignFineInv(raw, divisors)
hInv = max(1, ceil(raw));
for d = divisors(:).'
    if d <= 0, continue; end
    hInv = ceil(hInv / d) * d;
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


function s = gammaText(gamma)
if isempty(gamma)
    s = 'CIP gamma defaults from assembleCIP2D.';
else
    s = sprintf('CIP gamma=%g.', gamma);
end
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
