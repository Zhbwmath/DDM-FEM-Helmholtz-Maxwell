% VERIFY_LOD_PETERSEIM_FIG56_SQUARE  Square-domain surrogate for Peterseim Fig. 5-6.
%
% This reproduces the Figure 5-6 study pattern on Omega=(0,1)^2 with
% homogeneous impedance boundary conditions on the whole boundary. It is not
% the paper's triangle-scattering geometry.

fprintf('========== Peterseim Figure 5-6 Square-Domain LOD Verification ==========\n\n');

cfg = squareFig56Config();
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
outDir = fullfile(repoRoot, 'tasks', 'LOD', 'LOD_Helmholtz');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('Configuration: p=%d, h=1/%d, H^{-1}=[%s], k=[%s]\n', ...
    cfg.degree, cfg.hInv, num2str(cfg.hSweepInv), num2str(cfg.kSweep));
fprintf('Fine-scale rule: h <= C k_max^{-(2p+1)/(2p)}, k_max=%g, required h^{-1} >= %.3g\n', ...
    cfg.kMax, cfg.requiredHInv);
fprintf('Figure 6 surrogate: k=%g, ell=[%s]\n\n', cfg.kFixed, num2str(cfg.ellSweep));

fig5 = runFigure5Square(cfg);
fig6 = runFigure6Square(cfg);

fig5Path = fullfile(outDir, 'fig_peterseim_fig5_square.png');
fig6Path = fullfile(outDir, 'fig_peterseim_fig6_square.png');
matPath = fullfile(outDir, 'peterseim_fig56_square_results.mat');
mdPath = fullfile(outDir, 'peterseim_fig56_square_results.md');

plotFigure5Square(fig5, cfg, fig5Path);
plotFigure6Square(fig6, cfg, fig6Path);
save(matPath, 'cfg', 'fig5', 'fig6');
writeFigure56SquareReport(mdPath, cfg, fig5, fig6, fig5Path, fig6Path);

assert(all(isfinite([fig5.errLOD(:); fig5.errBest(:); fig5.errStab(:); fig5.errP1(:)])), ...
    'Figure 5 square surrogate has non-finite errors.');
assert(all(isfinite([fig6.errLOD(:); fig6.errBest(:); fig6.errStab(:); fig6.errP1(:)])), ...
    'Figure 6 square surrogate has non-finite errors.');

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', fig5Path, fig6Path, matPath, mdPath);
fprintf('\n========== Peterseim Figure 5-6 square-domain verification complete ==========\n');


function cfg = squareFig56Config()
cfg = struct();
cfg.hSweepInv = envVector('LOD_FIG56_HSWEEP_INV', [2, 4, 8]);
cfg.kSweep = envVector('LOD_FIG56_KSWEEP', [4, 8, 16]);
cfg.kFixed = envNumber('LOD_FIG56_KFIXED', 16);
cfg.ellSweep = envVector('LOD_FIG56_ELLS', [1, 2, 3]);
cfg.hSweepInvFig6 = envVector('LOD_FIG56_HSWEEP_INV_FIG6', [4, 8]);
cfg.degree = envNumber('LOD_FIG56_DEGREE', 1);
cfg.fineResolutionConstant = envNumber('LOD_FIG56_FINE_C', 1);
cfg.kMax = max([cfg.kSweep(:); cfg.kFixed]);
cfg.requiredHInv = cfg.fineResolutionConstant * cfg.kMax^((2 * cfg.degree + 1) / (2 * cfg.degree));
cfg.hInv = nestedFineHinv(envNumber('LOD_FIG56_HINV', NaN), cfg.requiredHInv, ...
    [cfg.hSweepInv, cfg.hSweepInvFig6]);
cfg.useParfor = logical(envNumber('LOD_FIG56_PARFOR', 0));
cfg.source = 1;
cfg.boundaryData = 0;
end


function fig5 = runFigure5Square(cfg)
nH = numel(cfg.hSweepInv);
nK = numel(cfg.kSweep);
fig5 = initResult(nH, nK);

for ik = 1:nK
    k = cfg.kSweep(ik);
    for iH = 1:nH
        Hinv = cfg.hSweepInv(iH);
        ell = max(1, round(log2(Hinv)));
        fprintf('Figure 5 surrogate: k=%g, H=1/%d, ell=%d ... ', k, Hinv, ell);
        fig5 = setResult(fig5, iH, ik, runSquareCase(k, Hinv, cfg.hInv, ell, cfg));
        fprintf('LOD %.3e, P1 %.3e\n', fig5.errLOD(iH, ik), fig5.errP1(iH, ik));
    end
end
fig5.Hinv = cfg.hSweepInv(:);
fig5.k = cfg.kSweep(:);
end


function fig6 = runFigure6Square(cfg)
nH = numel(cfg.hSweepInvFig6);
nL = numel(cfg.ellSweep);
fig6 = initResult(nH, nL);

for il = 1:nL
    ell = cfg.ellSweep(il);
    for iH = 1:nH
        Hinv = cfg.hSweepInvFig6(iH);
        fprintf('Figure 6 surrogate: k=%g, H=1/%d, ell=%d ... ', cfg.kFixed, Hinv, ell);
        fig6 = setResult(fig6, iH, il, runSquareCase(cfg.kFixed, Hinv, cfg.hInv, ell, cfg));
        fprintf('LOD %.3e, P1 %.3e\n', fig6.errLOD(iH, il), fig6.errP1(iH, il));
    end
end
fig6.Hinv = cfg.hSweepInvFig6(:);
fig6.ell = cfg.ellSweep(:);
fig6.kFixed = cfg.kFixed;
end


function result = initResult(nRow, nCol)
z = NaN(nRow, nCol);
result = struct('errLOD', z, 'errBest', z, 'errStab', z, 'errP1', z, ...
    'coarseDof', z, 'fineDof', z, 'ell', z, 'elapsed', z);
end


function result = setResult(result, i, j, one)
result.errLOD(i, j) = one.errLOD;
result.errBest(i, j) = one.errBest;
result.errStab(i, j) = one.errStab;
result.errP1(i, j) = one.errP1;
result.coarseDof(i, j) = one.coarseDof;
result.fineDof(i, j) = one.fineDof;
result.ell(i, j) = one.ell;
result.elapsed(i, j) = one.elapsed;
end


function one = runSquareCase(k, Hinv, hInv, ell, cfg)
timer = tic;
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1 / Hinv);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 1 / hInv);
opts = struct('oversampling', ell, 'solveCoarse', true, 'useParfor', cfg.useParfor);
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, nodeh, elemh, bdh, ...
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

[AH, bH] = assembleHelmholtz2D(nodeH, elemH, bdH, k, cfg.source, cfg.boundaryData, 1);
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
one.ell = ell;
one.elapsed = toc(timer);
end


function e = relEnergy(uf, uh, E, fineNorm)
d = uf - uh;
e = sqrt(max(0, real(d' * E * d))) / fineNorm;
end


function plotFigure5Square(fig5, cfg, outPath)
fig = figure('Name', 'Peterseim Fig5 Square Surrogate', 'Visible', 'off', ...
    'Position', [80, 80, 980, 760]);
x = 1 ./ fig5.Hinv;
titles = {'LOD Petrov-Galerkin', 'Best in LOD trial space', ...
    '$V_H$ trial / corrected test', 'Standard P1 FEM'};
data = {fig5.errLOD, fig5.errBest, fig5.errStab, fig5.errP1};
for p = 1:4
    subplot(2, 2, p);
    hold on;
    for ik = 1:numel(fig5.k)
        loglog(x, data{p}(:, ik), '-o', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('$\\kappa=%g$', fig5.k(ik)));
    end
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('$H$', 'Interpreter', 'latex');
    ylabel('relative $V$-error', 'Interpreter', 'latex');
    title(titles{p}, 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex');
end
sgtitle(sprintf('Square-domain surrogate of Peterseim Fig. 5 ($h=1/%d$, $\\ell=|\\log_2 H|$)', cfg.hInv), ...
    'Interpreter', 'latex');
print(fig, outPath, '-dpng', '-r200');
close(fig);
end


function plotFigure6Square(fig6, cfg, outPath)
fig = figure('Name', 'Peterseim Fig6 Square Surrogate', 'Visible', 'off', ...
    'Position', [80, 80, 1080, 360]);
x = 1 ./ fig6.Hinv;
titles = {'LOD Petrov-Galerkin', 'Best in LOD trial space', ...
    '$V_H$ trial / corrected test'};
data = {fig6.errLOD, fig6.errBest, fig6.errStab};
for p = 1:3
    subplot(1, 3, p);
    hold on;
    for il = 1:numel(fig6.ell)
        loglog(x, data{p}(:, il), '-o', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('$\\ell=%d$', fig6.ell(il)));
    end
    loglog(x, fig6.errP1(:, 1), 'k--s', 'LineWidth', 1.0, ...
        'DisplayName', 'P1 FEM');
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('$H$', 'Interpreter', 'latex');
    ylabel('relative $V$-error', 'Interpreter', 'latex');
    title(titles{p}, 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex');
end
sgtitle(sprintf('Square-domain surrogate of Peterseim Fig. 6 ($h=1/%d$, $\\kappa=%g$)', cfg.hInv, cfg.kFixed), ...
    'Interpreter', 'latex');
print(fig, outPath, '-dpng', '-r200');
close(fig);
end


function writeFigure56SquareReport(path, cfg, fig5, fig6, fig5Path, fig6Path)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Reproduction target: Peterseim Figures 5-6 square-domain surrogate.\n');
fprintf(fid, 'Created: 2026-05-26\n');
fprintf(fid, 'Updated: 2026-05-26\n');
fprintf(fid, 'Verification entry point: `verify/verify_lod_peterseim_fig56_square.m`\n');
fprintf(fid, 'Main utilities: `buildLODHelmholtz2D`, `weightedClementP1`, `assembleHelmholtz2D`\n\n');
fprintf(fid, '# Peterseim Figure 5-6 Square-Domain Surrogate\n\n');
fprintf(fid, 'This run follows the Figure 5-6 error-study pattern from Peterseim, but uses the normal square domain `Omega=(0,1)^2` with homogeneous impedance boundary conditions on the whole boundary. It does not include the paper''s triangular scatterer.\n\n');
fprintf(fid, 'Fine mesh rule: Lagrange degree `p=%d`, `k_max=%g`, require `h^{-1} >= %.6g` from `h = O(k^{-(2p+1)/(2p)})`; nested run uses `h=1/%d`.\n\n', ...
    cfg.degree, cfg.kMax, cfg.requiredHInv, cfg.hInv);
fprintf(fid, 'Default run: `h=1/%d`, Figure 5 `H^{-1}=[%s]`, `k=[%s]`; Figure 6 `k=%g`, `ell=[%s]`.\n\n', ...
    cfg.hInv, num2str(cfg.hSweepInv), num2str(cfg.kSweep), cfg.kFixed, num2str(cfg.ellSweep));
fprintf(fid, '![Figure 5 square surrogate](%s)\n\n', filepartsToMd(fig5Path));
fprintf(fid, '![Figure 6 square surrogate](%s)\n\n', filepartsToMd(fig6Path));

fprintf(fid, '## Figure 5 Data\n\n');
fprintf(fid, '| k | H^{-1} | ell | LOD | best LOD | stabilized | P1 | seconds |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ik = 1:numel(fig5.k)
    for iH = 1:numel(fig5.Hinv)
        fprintf(fid, '| %.8g | %d | %d | %.6e | %.6e | %.6e | %.6e | %.2f |\n', ...
            fig5.k(ik), fig5.Hinv(iH), fig5.ell(iH, ik), fig5.errLOD(iH, ik), ...
            fig5.errBest(iH, ik), fig5.errStab(iH, ik), fig5.errP1(iH, ik), ...
            fig5.elapsed(iH, ik));
    end
end

fprintf(fid, '\n## Figure 6 Data\n\n');
fprintf(fid, '| ell | H^{-1} | LOD | best LOD | stabilized | P1 | seconds |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|\n');
for il = 1:numel(fig6.ell)
    for iH = 1:numel(fig6.Hinv)
        fprintf(fid, '| %d | %d | %.6e | %.6e | %.6e | %.6e | %.2f |\n', ...
            fig6.ell(il), fig6.Hinv(iH), fig6.errLOD(iH, il), ...
            fig6.errBest(iH, il), fig6.errStab(iH, il), fig6.errP1(iH, il), ...
            fig6.elapsed(iH, il));
    end
end
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
        error('verify_lod_peterseim_fig56_square:env', ...
            'Environment variable %s must be numeric.', name);
    end
end
end


function val = envVector(name, defaultVal)
raw = getenv(name);
if isempty(raw)
    val = defaultVal;
else
    val = str2num(raw); %#ok<ST2NM>
    if isempty(val) || ~isnumeric(val)
        error('verify_lod_peterseim_fig56_square:env', ...
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
