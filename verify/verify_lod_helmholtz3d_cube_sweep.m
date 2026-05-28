% VERIFY_LOD_HELMHOLTZ3D_CUBE_SWEEP  Small cube-domain 3D Helmholtz LOD sweep.

fprintf('========== LOD Helmholtz 3D Cube Mini-Sweep ==========\n\n');

cfg = cube3DConfig();
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
outDir = fullfile(repoRoot, 'tasks', 'LOD', 'LOD_Helmholtz');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('Configuration: p=%d, h=1/%d, H^{-1}=[%s], k=[%s], ell=[%s]\n', ...
    cfg.degree, cfg.hInv, num2str(cfg.hSweepInv), num2str(cfg.kSweep), num2str(cfg.ellSweep));
fprintf('Fine-scale rule: h <= C k_max^{-(2p+1)/(2p)}, k_max=%g, required h^{-1} >= %.3g\n\n', ...
    cfg.kMax, cfg.requiredHInv);

results = runCube3DSweep(cfg);

figPath = fullfile(outDir, 'fig_lod_helmholtz3d_cube.png');
logEllFigPath = fullfile(outDir, 'fig_lod_helmholtz3d_logell.png');
matPath = fullfile(outDir, 'lod_helmholtz3d_cube_results.mat');
mdPath = fullfile(outDir, 'lod_helmholtz3d_cube_results.md');

plotCube3DSweep(results, cfg, figPath);
plotCube3DLogEllTrend(results, cfg, logEllFigPath);
save(matPath, 'cfg', 'results');
writeCube3DReport(mdPath, cfg, results, figPath, logEllFigPath);

assert(all(isfinite([results.errLOD(:); results.errBest(:); results.errStab(:); results.errP1(:)])), ...
    '3D cube sweep has non-finite errors.');

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', figPath, logEllFigPath, matPath, mdPath);
fprintf('\n========== LOD Helmholtz 3D cube mini-sweep complete ==========\n');


function cfg = cube3DConfig()
cfg = struct();
cfg.hSweepInv = envVector('LOD_3D_HSWEEP_INV', [2, 4, 8]);
cfg.kSweep = envVector('LOD_3D_KSWEEP', [1, 2]);
defaultEll = unique(max(1, ceil(log2(cfg.hSweepInv))));
cfg.ellSweep = envVector('LOD_3D_ELLS', defaultEll);
cfg.ellPolicy = envString('LOD_3D_ELL_POLICY', 'cross-log2-hinv');
cfg.degree = envNumber('LOD_3D_DEGREE', 1);
cfg.fineResolutionConstant = envNumber('LOD_3D_FINE_C', 1);
cfg.kMax = max(cfg.kSweep(:));
cfg.requiredHInv = 2 * cfg.fineResolutionConstant * cfg.kMax^((2 * cfg.degree + 1) / (2 * cfg.degree));
cfg.hInv = nestedFineHinv(envNumber('LOD_3D_HINV', 16), cfg.requiredHInv, cfg.hSweepInv);
cfg.useParfor = logical(envNumber('LOD_3D_PARFOR', 0));
cfg.solverMode = envString('LOD_3D_SOLVER', 'direct');
cfg.source = 1;
cfg.boundaryData = 0;
end


function results = runCube3DSweep(cfg)
nH = numel(cfg.hSweepInv);
nK = numel(cfg.kSweep);
nL = numel(cfg.ellSweep);
z = NaN(nH, nK, nL);
results = struct('errLOD', z, 'errBest', z, 'errStab', z, 'errP1', z, ...
    'coarseDof', z, 'fineDof', z, 'fullPatchCount', z, ...
    'meanPatchFraction', z, 'elapsed', z, ...
    'Hinv', cfg.hSweepInv(:), 'k', cfg.kSweep(:), 'ell', cfg.ellSweep(:));

for ik = 1:nK
    k = cfg.kSweep(ik);
    for iH = 1:nH
        Hinv = cfg.hSweepInv(iH);
        for il = 1:nL
            ell = cfg.ellSweep(il);
            fprintf('3D cube: k=%g, H=1/%d, ell=%d ... ', k, Hinv, ell);
            one = runCube3DCase(k, Hinv, cfg.hInv, ell, cfg);
            results.errLOD(iH, ik, il) = one.errLOD;
            results.errBest(iH, ik, il) = one.errBest;
            results.errStab(iH, ik, il) = one.errStab;
            results.errP1(iH, ik, il) = one.errP1;
            results.coarseDof(iH, ik, il) = one.coarseDof;
            results.fineDof(iH, ik, il) = one.fineDof;
            results.fullPatchCount(iH, ik, il) = one.fullPatchCount;
            results.meanPatchFraction(iH, ik, il) = one.meanPatchFraction;
            results.elapsed(iH, ik, il) = one.elapsed;
            fprintf('LOD %.3e, P1 %.3e, full patches %d/%d\n', ...
                one.errLOD, one.errP1, one.fullPatchCount, size(elemHForCount(Hinv), 1));
        end
    end
end
end


function one = runCube3DCase(k, Hinv, hInv, ell, cfg)
timer = tic;
[nodeH, elemH, bdH] = cubemesh([0, 1, 0, 1, 0, 1], 1 / Hinv);
[nodeh, elemh, bdh] = cubemesh([0, 1, 0, 1, 0, 1], 1 / hInv);
opts = struct('oversampling', ell, 'solveCoarse', true, ...
    'useParfor', cfg.useParfor, 'solverMode', cfg.solverMode);
lod = buildLODHelmholtz3D(nodeH, elemH, bdH, nodeh, elemh, bdh, ...
    k, cfg.source, cfg.boundaryData, opts);

A = lod.system.A;
b = lod.system.b;
E = lod.system.energy;
uf = lod.solution.fine;
fineNorm = sqrt(real(uf' * E * uf));
fineNorm = max(1, fineNorm);

P = lod.basis.coarse;
Psi = lod.basis.trial;
PsiStar = lod.basis.test;

[AH, bH] = assembleHelmholtz3D(nodeH, elemH, bdH, k, cfg.source, cfg.boundaryData, 1);
uP1 = P * (AH \ bH);

cBest = (Psi' * E * Psi) \ (Psi' * E * uf);
uBest = Psi * cBest;

cStab = (PsiStar' * A * P) \ (PsiStar' * b);
uStab = P * cStab;

one = struct();
one.errLOD = relEnergy(uf, lod.solution.uh, E, fineNorm);
one.errBest = relEnergy(uf, uBest, E, fineNorm);
one.errStab = relEnergy(uf, uStab, E, fineNorm);
one.errP1 = relEnergy(uf, uP1, E, fineNorm);
one.coarseDof = size(nodeH, 1);
one.fineDof = size(nodeh, 1);
patchElemCounts = cellfun(@numel, lod.patch.fineElemIds);
one.fullPatchCount = nnz(patchElemCounts == size(elemh, 1));
one.meanPatchFraction = mean(patchElemCounts) / size(elemh, 1);
one.elapsed = toc(timer);
end


function elemH = elemHForCount(Hinv)
[~, elemH, ~] = cubemesh([0, 1, 0, 1, 0, 1], 1 / Hinv);
end


function e = relEnergy(uf, uh, E, fineNorm)
d = uf - uh;
e = sqrt(max(0, real(d' * E * d))) / fineNorm;
end


function plotCube3DSweep(results, cfg, outPath)
fig = figure('Name', 'LOD Helmholtz 3D Cube Mini-Sweep', 'Visible', 'off', ...
    'Position', [80, 80, 960, 420]);
tiledlayout(1, numel(results.ell));
for il = 1:numel(results.ell)
    nexttile;
    hold on;
    for ik = 1:numel(results.k)
        loglog(1 ./ results.Hinv, results.errLOD(:, ik, il), '-o', ...
            'LineWidth', 1.2, 'DisplayName', sprintf('LOD k=%g', results.k(ik)));
        loglog(1 ./ results.Hinv, results.errP1(:, ik, il), '--s', ...
            'LineWidth', 1.0, 'DisplayName', sprintf('P1 k=%g', results.k(ik)));
    end
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('$H$', 'Interpreter', 'latex');
    ylabel('relative $V$-error', 'Interpreter', 'latex');
    title(sprintf('$\\ell=%d$', results.ell(il)), 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex');
end
sgtitle(sprintf('3D cube Helmholtz LOD mini-sweep ($h=1/%d$)', cfg.hInv), ...
    'Interpreter', 'latex');
print(fig, outPath, '-dpng', '-r200');
close(fig);
end


function plotCube3DLogEllTrend(results, cfg, outPath)
fig = figure('Name', 'LOD Helmholtz 3D Log Oversampling Trend', 'Visible', 'off', ...
    'Position', [80, 80, 980, 430]);
tiledlayout(1, numel(results.k));
H = 1 ./ results.Hinv(:);
ellForH = max(1, ceil(log2(results.Hinv(:))));
for ik = 1:numel(results.k)
    lodErr = NaN(size(results.Hinv));
    p1Err = NaN(size(results.Hinv));
    for iH = 1:numel(results.Hinv)
        il = find(results.ell == ellForH(iH), 1);
        if isempty(il)
            continue;
        end
        lodErr(iH) = results.errLOD(iH, ik, il);
        p1Err(iH) = results.errP1(iH, ik, il);
    end

    nexttile;
    loglog(H, lodErr, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', '$\mathrm{LOD},\ \ell=\lceil\log_2(H^{-1})\rceil$');
    hold on;
    loglog(H, p1Err, '--s', 'LineWidth', 1.2, 'MarkerSize', 6, ...
        'DisplayName', 'P1 FEM reference');
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('$H$ (decreases to the right)', 'Interpreter', 'latex');
    ylabel('relative $V$-error', 'Interpreter', 'latex');
    title(sprintf('k = %g', results.k(ik)));
    legend('Location', 'best', 'Interpreter', 'latex');
end
sgtitle(sprintf('3D cube Helmholtz LOD with logarithmic oversampling (h = 1/%d)', cfg.hInv));
print(fig, outPath, '-dpng', '-r200');
close(fig);
end


function writeCube3DReport(path, cfg, results, figPath, logEllFigPath)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: 3D cube Helmholtz LOD construction mini-sweep.\n');
fprintf(fid, 'Created: 2026-05-27\n');
fprintf(fid, 'Updated: 2026-05-27\n');
fprintf(fid, 'Verification entry point: `verify/verify_lod_helmholtz3d_cube_sweep.m`\n');
fprintf(fid, 'Main utilities: `buildLODHelmholtz3D`, `weightedClementP1`, `assembleHelmholtz3D`\n\n');
fprintf(fid, '# LOD Helmholtz 3D Cube Mini-Sweep\n\n');
fprintf(fid, 'This run checks the 3D tetrahedral LOD construction on `Omega=(0,1)^3` with homogeneous impedance boundary conditions. It is a small construction and behavior check, not a full Peterseim reproduction.\n\n');
fprintf(fid, 'Fine mesh rule: Lagrange degree `p=%d`, `k_max=%g`, require `h^{-1} >= %.6g` from `h = O(k^{-(2p+1)/(2p)})`; nested run uses `h=1/%d`.\n\n', ...
    cfg.degree, cfg.kMax, cfg.requiredHInv, cfg.hInv);
fprintf(fid, 'Default run: `h=1/%d`, `H^{-1}=[%s]`, `k=[%s]`, `ell=[%s]`.\n\n', ...
    cfg.hInv, num2str(cfg.hSweepInv), num2str(cfg.kSweep), num2str(cfg.ellSweep));
fprintf(fid, 'Corrector solver mode: `%s`. `direct` uses MATLAB sparse backslash for the local constrained saddle systems; `lu` uses explicit `decomposition(...,''lu'')` factor reuse for paired primal/adjoint solves.\n\n', cfg.solverMode);
fprintf(fid, 'Oversampling policy: `%s`. The default `ell` list is derived from `ceil(log2(H^{-1}))`; the full table is still combinatorial so the diagonal rows `ell=ceil(log2(H^{-1}))` can be compared against larger oversampling.\n\n', cfg.ellPolicy);
fprintf(fid, 'Patch coverage diagnostic: fixed `ell` does not imply comparable effective localization across different `H`. Use the full patch count and mean patch fraction columns to identify rows that are global-corrector cases or rows whose localized patches are much smaller than the coarser-mesh comparison.\n\n');
fprintf(fid, '## Log-Oversampling Trend\n\n');
fprintf(fid, 'The following figure compares the diagonal choice `ell=ceil(log2(H^{-1}))` against the P1 FEM reference as `H` decreases.\n\n');
fprintf(fid, '![3D cube logarithmic oversampling trend](%s)\n\n', filepartsToMd(logEllFigPath));

fprintf(fid, '## Full Combinatorial Sweep\n\n');
fprintf(fid, '![3D cube mini-sweep](%s)\n\n', filepartsToMd(figPath));

fprintf(fid, '## Data\n\n');
fprintf(fid, '| k | H^{-1} | ell | LOD | best LOD | stabilized | P1 | full patches | mean patch fraction | seconds |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ik = 1:numel(results.k)
    for iH = 1:numel(results.Hinv)
        for il = 1:numel(results.ell)
            fprintf(fid, '| %.8g | %d | %d | %.6e | %.6e | %.6e | %.6e | %d/%d | %.3f | %.2f |\n', ...
                results.k(ik), results.Hinv(iH), results.ell(il), ...
                results.errLOD(iH, ik, il), results.errBest(iH, ik, il), ...
                results.errStab(iH, ik, il), results.errP1(iH, ik, il), ...
                results.fullPatchCount(iH, ik, il), coarseElementCount(results.Hinv(iH)), ...
                results.meanPatchFraction(iH, ik, il), ...
                results.elapsed(iH, ik, il));
        end
    end
end
end


function n = coarseElementCount(Hinv)
[~, elem, ~] = cubemesh([0, 1, 0, 1, 0, 1], 1 / Hinv);
n = size(elem, 1);
end


function s = filepartsToMd(path)
[~, name, ext] = fileparts(path);
s = [name, ext];
end


function val = envNumber(name, defaultVal)
raw = getenv(name);
if isempty(raw)
    val = defaultVal;
else
    val = str2double(raw);
    if ~isfinite(val)
        error('verify_lod_helmholtz3d_cube_sweep:env', ...
            'Environment variable %s must be numeric.', name);
    end
end
end


function val = envString(name, defaultVal)
raw = getenv(name);
if isempty(raw)
    val = defaultVal;
else
    val = raw;
end
end


function val = envVector(name, defaultVal)
raw = getenv(name);
if isempty(raw)
    val = defaultVal;
else
    val = str2num(raw); %#ok<ST2NM>
    if isempty(val) || ~isnumeric(val)
        error('verify_lod_helmholtz3d_cube_sweep:env', ...
            'Environment variable %s must be a numeric vector.', name);
    end
end
val = val(:).';
end


function hInv = nestedFineHinv(requested, required, coarseInvs)
base = 1;
coarseInvs = unique(round(coarseInvs(:).'));
for j = 1:numel(coarseInvs)
    base = lcm(base, coarseInvs(j));
end
if isfinite(requested)
    target = max(requested, required);
else
    target = required;
end
hInv = base * ceil(target / base);
end
