% VERIFY_DDM3LVL_LOD_HELMHOLTZ_EXPERIMENTS  Gated three-level LOD-DDM sweep.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
outDir = fullfile(repoRoot, 'tasks', 'DDM3lvl_LOD_Helmholtz');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg = config();
cases = experimentCases(cfg);
results = repmat(emptyResult(), numel(cases), 1);
csvPath = fullfile(outDir, ['three_level_results' cfg.outputSuffix '.csv']);
mdPath = fullfile(outDir, ['three_level_results' cfg.outputSuffix '.md']);

fprintf('========== Three-Level LOD-DDM Helmholtz Experiments ==========\n');
fprintf('Cases: %d, run enabled: %d, max run DOF: %d\n\n', ...
    numel(cases), cfg.runEnabled, cfg.maxRunDof);

for i = 1:numel(cases)
    c = finalizeCase(cases(i));
    r = copyCaseFields(emptyResult(), c);
    est = estimateCase(c, cfg);
    r.ndof = est.N;
    r.ncoarse = est.Nc;
    r.estimateGB = est.totalGB;

    fprintf('[%02d/%02d] k=%g eps=%s H0=1/%d m=%d: N=%d estimate %.2f GB ... ', ...
        i, numel(cases), c.k, c.epsilonLabel, c.H0Inv, c.m, est.N, est.totalGB);

    if est.totalGB > cfg.memoryLimitGB
        r.status = 'blocked_memory_gt_limit';
        r.notes = sprintf('Estimate %.2f GB exceeds %.2f GB.', est.totalGB, cfg.memoryLimitGB);
        fprintf('BLOCKED\n');
    elseif ~cfg.runEnabled
        r.status = 'estimated_only';
        r.notes = 'Set DDM3LVL_RUN=1 to execute memory-permitted rows.';
        fprintf('ESTIMATED\n');
    elseif est.N > cfg.maxRunDof
        r.status = 'queued_runtime_cap';
        r.notes = sprintf('Below memory gate but above max run DOF %d.', cfg.maxRunDof);
        fprintf('QUEUED\n');
    else
        try
            out = runOneCase(c, cfg);
            r = mergeStruct(r, out);
            r.status = 'ran';
            fprintf('RAN (exact %d, three %d, coarse %d)\n', ...
                r.exactIter, r.threeIter, r.coarseIter);
        catch ME
            r.status = 'failed';
            r.notes = sprintf('%s: %s', ME.identifier, ME.message);
            fprintf('FAILED: %s\n', ME.message);
        end
    end
    results(i) = r;
    writeOutputs(csvPath, mdPath, results(1:i), cfg);
end

writeOutputs(csvPath, mdPath, results, cfg);
writePlots(outDir, results, cfg);
fprintf('\nCSV: %s\n', csvPath);
fprintf('Markdown: %s\n', mdPath);
fprintf('========== Three-Level LOD-DDM experiment pass complete ==========\n');


function cfg = config()
cfg.kValues = envVector('DDM3LVL_K_VALUES', [16, 32, 64]);
cfg.H0InvValues = envVector('DDM3LVL_H0_INV_VALUES', [1, 2, 4, 8]);
cfg.mValues = envVector('DDM3LVL_M_VALUES', [1, 2, 3, 4]);
cfg.epsilonLabels = {'zero', 'k'};
cfg.runEnabled = logical(envNumber('DDM3LVL_RUN', 0));
cfg.maxRunDof = envNumber('DDM3LVL_MAX_RUN_DOF', 1500);
cfg.memoryLimitGB = envNumber('DDM3LVL_MEMORY_LIMIT_GB', 200);
cfg.tol = envNumber('DDM3LVL_TOL', 1e-6);
cfg.maxit = envNumber('DDM3LVL_MAXIT', 80);
cfg.smax = envNumber('DDM3LVL_SMAX', 12);
cfg.Ch = envNumber('DDM3LVL_CH', 1);
suffix = envString('DDM3LVL_OUTPUT_SUFFIX', '');
if ~isempty(suffix) && suffix(1) ~= '_'
    suffix = ['_' suffix];
end
cfg.outputSuffix = suffix;
end


function cases = experimentCases(cfg)
cases = repmat(baseCase(), 0, 1);
for k = cfg.kValues
    for H0Inv = cfg.H0InvValues
        for m = cfg.mValues
            for e = 1:numel(cfg.epsilonLabels)
                c = baseCase();
                c.k = k;
                c.H0Inv = H0Inv;
                c.m = m;
                c.epsilonLabel = cfg.epsilonLabels{e};
                cases(end+1) = c; %#ok<AGROW>
            end
        end
    end
end
end


function c = baseCase()
c = struct('k', NaN, 'epsilonLabel', '', 'hInv', NaN, 'HInv', NaN, ...
    'H0Inv', NaN, 'm', NaN);
end


function c = finalizeCase(c)
c.hInv = alignFineInv(ceil(c.k^(3/2)), [c.k, c.H0Inv]);
c.HInv = compatibleDivisor(c.hInv, max(2, c.k / 2));
c.H0Inv = compatibleDivisor(c.hInv, c.H0Inv);
end


function est = estimateCase(c, cfg)
N = (c.hInv + 1)^2;
NT = 2 * c.hInv^2;
Nc = (c.HInv + 1)^2;
patchWidth = min(1, (2 * c.m + 1) / max(1, c.HInv));
nPatch = max(4, (ceil(patchWidth * c.hInv) + 1)^2);
globalBytes = 112 * N + 48 * NT;
gmresBytes = 16 * N * (cfg.maxit + 3);
lodBytes = 2 * 16 * Nc * min(N, nPatch);
coarseBytes = 112 * Nc + 16 * Nc * Nc;
est = struct('N', N, 'NT', NT, 'Nc', Nc, ...
    'totalGB', (globalBytes + gmresBytes + lodBytes + coarseBytes) / 2^30);
end


function out = runOneCase(c, cfg)
h = 1 / c.hInv;
H = 1 / c.HInv;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], h);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], H);
epsilon = epsilonValue(c.epsilonLabel);
pde = helmholtzPDE(c.k, 'epsilon', epsilon, 'eta', 'sqrt');
partsFine = coarseHatPartition2D(node, elem, bdFlag, 1 / max(1, c.H0Inv));
lodOpts = struct('oversampling', c.m, 'solveCoarse', false, 'solverMode', 'direct');
preOpts = struct('variant', 'impedance', 'coarseType', 'lod', ...
    'lodOptions', lodOpts, 'solverMode', 'direct', 'adjointType', 'energy');
preExact = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, pde, ...
    partsFine, nodeH, elemH, bdH, preOpts);

grid = max(1, [c.H0Inv, c.H0Inv]);
coarseOpts = struct('subdomainGrid', grid, 'overlap', 1 / max(1, c.H0Inv), ...
    'greaterOverlap', 1 / max(1, c.H0Inv), 'smax', cfg.smax, ...
    'compareLocalBasis', true);
coarse = buildLODCoarseSchwarzHelmholtz2D(preExact, nodeH, elemH, bdH, coarseOpts);
threeOpts = struct('fineSpace', preExact.fineSpace, ...
    'coarseSpace', coarse.coarseSpace, 'localSolver', preExact.localSolver, ...
    'variant', 'impedance', 'solverMode', 'direct', 'adjointType', 'energy');
preThree = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, pde, ...
    partsFine, nodeH, elemH, bdH, threeOpts);

b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, c.k, 1);
[~, flagExact, relExact, iterExact, resExact] = gmres(preExact.A, b, [], cfg.tol, cfg.maxit, @preExact.applyResidual);
[~, flagThree, relThree, iterThree, resThree] = gmres(preThree.A, b, [], cfg.tol, cfg.maxit, @preThree.applyResidual);
r0 = preExact.coarseSpace.test' * b;
[~, flagCoarse, relCoarse, iterCoarse, resCoarse] = gmres(coarse.A0, r0, [], cfg.tol, cfg.maxit, @coarse.applyM0inv);

out = struct();
out.exactIter = gmresIterationCount(iterExact);
out.threeIter = gmresIterationCount(iterThree);
out.coarseIter = gmresIterationCount(iterCoarse);
out.exactFlag = flagExact;
out.threeFlag = flagThree;
out.coarseFlag = flagCoarse;
out.exactRelres = relExact;
out.threeRelres = relThree;
out.coarseRelres = relCoarse;
out.sContract = coarse.diagnostics.sContract;
out.sFovPositive = coarse.diagnostics.sFovPositive;
out.maxNormEPower = max(coarse.diagnostics.normEPower);
out.maxAlpha = max(coarse.diagnostics.alpha);
out.basisContainedTrial = coarse.basisComparison.maxContainedTrialRelEnergy;
out.basisContainedTest = coarse.basisComparison.maxContainedTestRelEnergy;
out.basisKernelTrial = coarse.basisComparison.maxKernelTrial;
out.basisKernelTest = coarse.basisComparison.maxKernelTest;
out.exactResvecLength = numel(resExact);
out.threeResvecLength = numel(resThree);
out.coarseResvecLength = numel(resCoarse);
out.notes = '';
end


function eps = epsilonValue(label)
switch lower(label)
    case 'zero'
        eps = 0;
    case 'k'
        eps = 'k';
    otherwise
        error('verify_ddm3lvl:epsilon', 'Unknown epsilon label "%s".', label);
end
end


function r = emptyResult()
r = struct('k', NaN, 'epsilonLabel', '', 'hInv', NaN, 'HInv', NaN, ...
    'H0Inv', NaN, 'm', NaN, 'ndof', NaN, 'ncoarse', NaN, 'estimateGB', NaN, ...
    'exactIter', NaN, 'threeIter', NaN, 'coarseIter', NaN, ...
    'exactFlag', NaN, 'threeFlag', NaN, 'coarseFlag', NaN, ...
    'exactRelres', NaN, 'threeRelres', NaN, 'coarseRelres', NaN, ...
    'sContract', NaN, 'sFovPositive', NaN, 'maxNormEPower', NaN, ...
    'maxAlpha', NaN, 'basisContainedTrial', NaN, 'basisContainedTest', NaN, ...
    'basisKernelTrial', NaN, 'basisKernelTest', NaN, ...
    'exactResvecLength', NaN, 'threeResvecLength', NaN, ...
    'coarseResvecLength', NaN, 'status', '', 'notes', '');
end


function r = copyCaseFields(r, c)
r.k = c.k;
r.epsilonLabel = c.epsilonLabel;
r.hInv = c.hInv;
r.HInv = c.HInv;
r.H0Inv = c.H0Inv;
r.m = c.m;
end


function r = mergeStruct(r, out)
names = fieldnames(out);
for i = 1:numel(names)
    r.(names{i}) = out.(names{i});
end
end


function writeOutputs(csvPath, mdPath, results, cfg)
writeCsv(csvPath, results);
writeMarkdown(mdPath, results, cfg);
end


function writeCsv(path, results)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
names = fieldnames(results);
fprintf(fid, '%s\n', strjoin(names.', ','));
for i = 1:numel(results)
    vals = cell(1, numel(names));
    for j = 1:numel(names)
        v = results(i).(names{j});
        if isnumeric(v)
            vals{j} = numericCsvValue(v);
        else
            vals{j} = csvEscape(v);
        end
    end
    fprintf(fid, '%s\n', strjoin(vals, ','));
end
end


function writeMarkdown(path, results, cfg)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.\n');
fprintf(fid, 'Created: 2026-06-18\n');
fprintf(fid, 'Updated: %s\n', char(datetime('today', 'Format', 'yyyy-MM-dd')));
fprintf(fid, 'Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`\n');
fprintf(fid, 'Main utilities: `buildLODCoarseSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`\n\n');
fprintf(fid, '# Three-Level LOD-DDM Experiment Results\n\n');
fprintf(fid, 'Run enabled: `%d`. Memory limit: %.1f GB. Max run DOF: %d.\n\n', ...
    cfg.runEnabled, cfg.memoryLimitGB, cfg.maxRunDof);
fprintf(fid, '| k | eps | h | H | H0 | m | N | estimate GB | exact it | three it | coarse it | s0 | alpha s | basis contained trial | basis contained test | status | notes |\n');
fprintf(fid, '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '| %.0f | %s | 1/%d | 1/%d | 1/%d | %d | %.0f | %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
        r.k, r.epsilonLabel, r.hInv, r.HInv, r.H0Inv, r.m, r.ndof, r.estimateGB, ...
        mdNumber(r.exactIter), mdNumber(r.threeIter), mdNumber(r.coarseIter), ...
        mdNumber(r.sContract), mdNumber(r.sFovPositive), ...
        mdNumber(r.basisContainedTrial), mdNumber(r.basisContainedTest), ...
        r.status, mdEscape(r.notes));
end
end


function writePlots(outDir, results, cfg)
ran = strcmp({results.status}, 'ran');
if ~any(ran)
    return;
end
r = results(ran);
fig = figure('Visible', 'off');
plot([r.k], [r.exactIter], 'o-', 'LineWidth', 1.2); hold on;
plot([r.k], [r.threeIter], 's-', 'LineWidth', 1.2);
xlabel('$k$', 'Interpreter', 'latex');
ylabel('GMRES iterations', 'Interpreter', 'latex');
legend({'exact two-level', 'three-level'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveas(fig, fullfile(outDir, ['fig_three_level_gmres_iterations' cfg.outputSuffix '.png']));
close(fig);

fig = figure('Visible', 'off');
plot([r.k], [r.sContract], 'o-', 'LineWidth', 1.2);
xlabel('$k$', 'Interpreter', 'latex');
ylabel('$s_0$', 'Interpreter', 'latex');
grid on;
saveas(fig, fullfile(outDir, ['fig_three_level_s0_vs_k' cfg.outputSuffix '.png']));
close(fig);
end


function n = alignFineInv(raw, targetInvs)
targetInvs = round(targetInvs(isfinite(targetInvs) & targetInvs > 0));
base = 1;
for i = 1:numel(targetInvs)
    base = lcm(base, max(1, targetInvs(i)));
end
n = ceil(raw / base) * base;
end


function d = compatibleDivisor(n, target)
target = max(1, round(target));
divs = divisorsInt(n);
[~, idx] = min(abs(divs - target));
d = divs(idx);
end


function d = divisorsInt(n)
d = [];
for q = 1:floor(sqrt(n))
    if mod(n, q) == 0
        d(end+1) = q; %#ok<AGROW>
        if q ~= n / q
            d(end+1) = n / q; %#ok<AGROW>
        end
    end
end
d = sort(d);
end


function n = gmresIterationCount(iter)
if isempty(iter)
    n = NaN;
elseif numel(iter) == 1
    n = iter;
else
    n = iter(2);
end
end


function v = envNumber(name, defaultValue)
s = getenv(name);
if isempty(s)
    v = defaultValue;
else
    v = str2double(s);
    if isnan(v), v = defaultValue; end
end
end


function v = envString(name, defaultValue)
s = getenv(name);
if isempty(s)
    v = defaultValue;
else
    v = s;
end
end


function v = envVector(name, defaultValue)
s = getenv(name);
if isempty(s)
    v = defaultValue;
else
    v = str2num(s); %#ok<ST2NM>
    if isempty(v), v = defaultValue; end
end
end


function s = csvEscape(v)
s = char(v);
if contains(s, ',') || contains(s, '"') || contains(s, newline)
    s = ['"', strrep(s, '"', '""'), '"'];
end
end


function s = numericCsvValue(v)
v = full(v);
if isempty(v)
    s = '';
elseif isscalar(v)
    s = sprintf('%.16g', v);
else
    parts = arrayfun(@(x) sprintf('%.16g', x), v(:).', 'UniformOutput', false);
    s = csvEscape(strjoin(parts, ' '));
end
end


function s = mdNumber(x)
x = full(x);
if ~isscalar(x)
    x = x(1);
end
if isnan(x)
    s = '-';
elseif abs(x) >= 1e4 || (abs(x) > 0 && abs(x) < 1e-3)
    s = sprintf('%.3e', x);
else
    s = sprintf('%.6g', x);
end
end


function s = mdEscape(v)
s = char(v);
s = strrep(s, '|', '\|');
s = strrep(s, newline, ' ');
end
