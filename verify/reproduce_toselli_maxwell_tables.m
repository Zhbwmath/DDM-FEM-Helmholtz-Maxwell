function results = reproduce_toselli_maxwell_tables(mode, opts)
% REPRODUCE_TOSELLI_MAXWELL_TABLES  Reproduce Toselli Maxwell ASM table cells.

if nargin < 1 || isempty(mode), mode = 'smoke'; end
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'tol'), opts.tol = 1e-6; end
if ~isfield(opts, 'maxit'), opts.maxit = 300; end
if ~isfield(opts, 'allowLarge'), opts.allowLarge = false; end
if ~isfield(opts, 'maxDefaultDofs'), opts.maxDefaultDofs = 50000; end
if ~isfield(opts, 'memoryLimitGB'), opts.memoryLimitGB = 200; end
if ~isfield(opts, 'writeReport'), opts.writeReport = true; end

cases = toselliCases(mode);
results = repmat(emptyResult(), numel(cases), 1);

fprintf('========== Toselli Maxwell Schwarz Reproduction (%s) ==========\n\n', mode);
for k = 1:numel(cases)
    c = cases(k);
    est = estimateToselliCase(c);
    results(k) = emptyResult(c, est);
    fprintf('Table %d n=%d m=%d H/delta=%.4g eta1=%.3g %s ... ', ...
        c.table, c.n, c.m, c.Hdelta, c.eta1, c.level);

    if est.deltaLayers ~= round(est.deltaLayers)
        results(k).status = "skipped";
        results(k).notes = "non-integer delta/h";
        fprintf('SKIPPED (%s)\n', results(k).notes);
        continue;
    end
    if est.memoryGB > opts.memoryLimitGB
        results(k).status = "blocked";
        results(k).notes = sprintf('estimated memory %.1f GB exceeds %.1f GB gate', ...
            est.memoryGB, opts.memoryLimitGB);
        fprintf('BLOCKED (%s)\n', results(k).notes);
        continue;
    end
    if ~opts.allowLarge && est.freeDofs > opts.maxDefaultDofs
        results(k).status = "blocked";
        results(k).notes = sprintf('free edge DOFs %d exceed default run limit %d', ...
            est.freeDofs, opts.maxDefaultDofs);
        fprintf('BLOCKED (%s)\n', results(k).notes);
        continue;
    end

    try
        r = runToselliCase(c, opts);
        results(k) = mergeResult(c, est, r);
        fprintf('DONE kappa %.2f (%d it), paper %.2f (%d)\n', ...
            r.condest, r.iter, c.paperKappa, c.paperIter);
    catch err
        results(k).status = "failed";
        results(k).notes = string(err.message);
        fprintf('FAILED (%s)\n', err.message);
    end
end

if opts.writeReport
    writeToselliReport(results, mode);
end
fprintf('\n========== Done ==========\n');
end


function r = runToselliCase(c, opts)
h = 1 / c.n;
H = 1 / c.m;
delta = H / c.Hdelta;

[node, elem, bdFlag] = cubemesh([0, 1, 0, 1, 0, 1], h);
[A, ~] = assembleMaxwell3D(node, elem, c.eta1, c.eta2);
[freeEdges, ~] = nedelecFreeEdges3D(elem, bdFlag);
A_ff = A(freeEdges, freeEdges);

parts = partitionMesh3D(node, elem, bdFlag, [c.m, c.m, c.m], 'overlap', delta);
edgeParts = nedelecSubdomainEdges3D(elem, bdFlag, parts);
[Mfun, info] = nedelecAdditiveSchwarz3D(A_ff, edgeParts, freeEdges);

if c.level == "two"
    [coarseNode, coarseElem, coarseBd] = cubemesh([0, 1, 0, 1, 0, 1], H);
    [Pfull, ~] = prolongateNestedNed1(coarseNode, coarseElem, node, elem);
    [coarseFree, ~] = nedelecFreeEdges3D(coarseElem, coarseBd);
    P_H = Pfull(freeEdges, coarseFree);
    [Mfun, info] = nedelecTwoLevelASM3D(A_ff, edgeParts, freeEdges, P_H);
end

b = deterministicRhs(size(A_ff, 1));
[~, stats] = pcgLanczosCondition(A_ff, b, opts.tol, opts.maxit, Mfun);

r = struct();
r.condest = stats.condest;
r.iter = stats.iter;
r.relres = stats.relres;
r.flag = stats.flag;
r.solverInfo = info;
end


function b = deterministicRhs(n)
j = (1:n).';
b = sin(0.37 * j) + cos(0.11 * j);
end


function est = estimateToselliCase(c)
H = 1 / c.m;
delta = H / c.Hdelta;
h = 1 / c.n;
freeDofs = round(6.82 * c.n^3);
nSub = c.m^3;
overlapScale = min((1 / c.m + 2 * delta)^3 / max((1 / c.m)^3, eps), c.m^3);
localDofs = max(1, freeDofs / nSub * overlapScale);
globalGB = 16 * 40 * freeDofs / 1e9;
localGB = nSub * 16 * 40 * localDofs * 6 / 1e9;
coarseGB = 0;
if c.level == "two"
    coarseDofs = round(6.82 * c.m^3);
    coarseGB = 16 * coarseDofs^2 / 1e9;
end
est = struct();
est.delta = delta;
est.deltaLayers = delta / h;
est.freeDofs = freeDofs;
est.memoryGB = globalGB + localGB + coarseGB;
end


function cases = toselliCases(mode)
switch lower(mode)
    case 'smoke'
        cases = [
            makeCase(1, "one", 8, 2, 4, 1, 1, 14.10, 23)
            makeCase(2, "two", 8, 2, 4, 1, 1, 8.94, 19)
            makeCase(3, "one", 16, 2, 4, 1, 1, 13.32, 21)
            makeCase(4, "two", 16, 2, 4, 1, 1, 8.49, 19)
        ];
    case 'tables34'
        cases = table34Cases();
    case 'full'
        cases = [table12Cases(); table34Cases()];
    otherwise
        error('reproduce_toselli_maxwell_tables:badMode', ...
            'Unknown mode %s. Use smoke, tables34, or full.', mode);
end
end


function cases = table12Cases()
rows = [8 2; 16 2; 16 4; 16 8; 24 3; 24 6; 24 12; 32 4; 32 8; 40 5; 40 10; 48 6; 48 12];
hd = [8 4 2 4/3];
t1k = [
    NaN 14.10 8.95 8.94
    28.05 13.32 NaN NaN
    NaN 33.51 15.13 27.30
    NaN NaN 57.31 NaN
    50.94 22.62 NaN NaN
    NaN 67.39 22.99 21.25
    NaN NaN 73.69 NaN
    80.86 34.82 NaN NaN
    NaN 117.18 37.33 28.85
    115.93 50.41 NaN NaN
    NaN 177.46 55.35 44.52
    143.05 66.60 NaN NaN
    NaN 243.85 74.57 60.03];
t1i = [
    NaN 23 20 15
    25 21 NaN NaN
    NaN 29 23 28
    NaN NaN 35 NaN
    30 24 NaN NaN
    NaN 35 26 27
    NaN NaN 38 NaN
    36 25 NaN NaN
    NaN 43 29 28
    41 29 NaN NaN
    NaN 46 31 29
    43 32 NaN NaN
    NaN 55 35 33];
t2k = [
    NaN 8.94 8.98 9.05
    14.05 8.49 NaN NaN
    NaN 8.61 9.86 27.54
    NaN NaN 10.59 NaN
    13.86 8.45 NaN NaN
    NaN 8.43 9.05 21.93
    NaN NaN 9.39 NaN
    13.02 8.30 NaN NaN
    NaN 8.37 8.78 21.39
    13.12 8.29 NaN NaN
    NaN 8.29 8.68 22.28
    12.91 8.36 NaN NaN
    NaN 8.32 8.64 22.93];
t2i = [
    NaN 19 20 15
    21 19 NaN NaN
    NaN 19 21 28
    NaN NaN 19 NaN
    21 19 NaN NaN
    NaN 19 19 27
    NaN NaN 18 NaN
    20 18 NaN NaN
    NaN 19 19 25
    20 18 NaN NaN
    NaN 18 19 25
    20 18 NaN NaN
    NaN 18 18 24];

cases = repmat(makeCase(1, "one", 8, 2, 4, 1, 1, NaN, NaN), 2*numel(rows)*numel(hd), 1);
idx = 0;
for r = 1:size(rows, 1)
    for q = 1:numel(hd)
        if ~isnan(t1k(r, q))
            idx = idx + 1;
            cases(idx) = makeCase(1, "one", rows(r,1), rows(r,2), hd(q), 1, 1, t1k(r,q), t1i(r,q));
        end
        if ~isnan(t2k(r, q))
            idx = idx + 1;
            cases(idx) = makeCase(2, "two", rows(r,1), rows(r,2), hd(q), 1, 1, t2k(r,q), t2i(r,q));
        end
    end
end
cases = cases(1:idx);
end


function cases = table34Cases()
eta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10 100 1000 10000 1e5 1e6];
defs = [
    table34Def(3, "one", 2, 8, ...
        [29.44 29.44 29.43 29.43 29.40 29.20 28.05 22.34 11.51 8.11 8.53 8.95 9.13], ...
        [26 26 26 26 26 26 25 22 17 14 13 14 14])
    table34Def(3, "one", 2, 4, ...
        [13.58 13.58 13.58 13.58 13.58 13.52 13.32 11.80 8.04 8.00 8.00 8.02 8.02], ...
        [22 22 22 22 22 22 21 18 15 12 10 9 9])
    table34Def(3, "one", 4, 4, ...
        [34.50 34.50 34.50 34.49 34.49 34.41 33.50 25.81 11.54 8.14 8.70 9.17 9.24], ...
        [31 31 31 31 31 31 29 25 18 14 14 15 15])
    table34Def(3, "one", 4, 2, ...
        [14.50 14.50 14.50 14.51 14.62 16.04 15.13 10.19 7.99 7.99 8.01 8.05 8.06], ...
        [25 25 25 25 25 25 23 20 17 14 12 10 10])
    table34Def(3, "one", 4, 4/3, ...
        [14.36 14.36 14.36 14.37 15.11 28.18 27.30 24.98 26.71 26.96 26.98 26.98 26.99], ...
        [24 24 24 24 25 28 28 29 26 25 23 22 22])
    table34Def(4, "two", 2, 8, ...
        [13.30 13.30 13.30 13.30 13.30 13.29 14.05 13.48 9.94 8.11 8.55 8.92 9.00], ...
        [21 21 21 21 21 21 21 19 16 15 15 14 14])
    table34Def(4, "two", 2, 4, ...
        [8.99 8.99 8.99 8.99 8.98 8.89 8.49 8.43 8.04 8.07 8.08 8.09 8.09], ...
        [20 20 20 20 20 20 19 17 15 13 13 12 12])
    table34Def(4, "two", 4, 4, ...
        [8.94 8.94 8.94 8.94 8.93 8.86 8.61 9.03 8.59 8.30 8.91 9.40 9.47], ...
        [20 20 20 20 20 20 19 18 18 17 17 18 18])
    table34Def(4, "two", 4, 2, ...
        [13.44 13.44 13.44 13.43 13.35 12.64 9.86 8.67 8.74 8.82 8.88 8.92 8.93], ...
        [24 24 24 24 24 23 21 19 18 17 15 14 14])
    table34Def(4, "two", 4, 4/3, ...
        [14.31 14.31 14.31 14.32 19.51 27.83 27.54 24.99 26.67 26.97 27.03 27.04 27.04], ...
        [24 24 24 24 26 27 28 29 26 26 23 24 24])
];

cases = repmat(makeCase(3, "one", 16, 2, 4, 1, 1, NaN, NaN), numel(defs)*numel(eta), 1);
idx = 0;
for d = 1:numel(defs)
    for j = 1:numel(eta)
        idx = idx + 1;
        cases(idx) = makeCase(defs(d).table, defs(d).level, 16, defs(d).m, ...
            defs(d).Hdelta, eta(j), 1, defs(d).kappa(j), defs(d).iter(j));
    end
end
end


function d = table34Def(tableId, level, m, Hdelta, kappa, iter)
d = struct('table', tableId, 'level', level, 'm', m, 'Hdelta', Hdelta, ...
    'kappa', kappa, 'iter', iter);
end


function c = makeCase(tableId, level, n, m, Hdelta, eta1, eta2, paperKappa, paperIter)
c = struct();
c.table = tableId;
c.level = level;
c.n = n;
c.m = m;
c.Hdelta = Hdelta;
c.eta1 = eta1;
c.eta2 = eta2;
c.paperKappa = paperKappa;
c.paperIter = paperIter;
end


function r = emptyResult(c, est)
if nargin < 1
    c = makeCase(NaN, "one", NaN, NaN, NaN, NaN, NaN, NaN, NaN);
end
if nargin < 2
    est = struct('delta', NaN, 'deltaLayers', NaN, 'freeDofs', NaN, 'memoryGB', NaN);
end
r = c;
r.delta = est.delta;
r.deltaLayers = est.deltaLayers;
r.freeDofsEstimate = est.freeDofs;
r.memoryGBEstimate = est.memoryGB;
r.repoKappa = NaN;
r.repoIter = NaN;
r.relres = NaN;
r.flag = NaN;
r.absKappaDiff = NaN;
r.relKappaDiff = NaN;
r.iterDiff = NaN;
r.status = "not-run";
r.notes = "";
end


function out = mergeResult(c, est, run)
out = emptyResult(c, est);
out.repoKappa = run.condest;
out.repoIter = run.iter;
out.relres = run.relres;
out.flag = run.flag;
out.absKappaDiff = abs(run.condest - c.paperKappa);
out.relKappaDiff = out.absKappaDiff / max(abs(c.paperKappa), eps);
out.iterDiff = run.iter - c.paperIter;
out.status = "done";
if run.flag ~= 0
    out.status = "failed";
    out.notes = sprintf('PCG flag %d', run.flag);
elseif isfinite(c.paperKappa)
    if out.relKappaDiff < 0.1
        out.notes = "matched";
    elseif out.relKappaDiff < 0.35
        out.notes = "close";
    else
        out.notes = "different";
    end
else
    out.notes = "paper target not encoded";
end
end


function writeToselliReport(results, mode)
outDir = fullfile('docs', 'Tos00_maxwell_schwarz');
if ~exist(outDir, 'dir'), mkdir(outDir); end
outPath = fullfile(outDir, sprintf('toselli_maxwell_%s_results.md', mode));
fid = fopen(outPath, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Reproduction target: Tables 1-4 in Toselli, *Overlapping Schwarz methods for Maxwell''s equations in three dimensions*.\n');
fprintf(fid, 'Created: 2026-05-26\n');
fprintf(fid, 'Updated: 2026-05-26\n');
fprintf(fid, 'Verification entry point: `verify/verify_toselli_maxwell_schwarz.m`; reproduction driver: `verify/reproduce_toselli_maxwell_tables.m`.\n');
fprintf(fid, 'Main utilities: `assembleMaxwell3D`, `nedelecAdditiveSchwarz3D`, `nedelecTwoLevelASM3D`, `prolongateNestedNed1`, `pcgLanczosCondition`.\n\n');
fprintf(fid, 'Mode: `%s`.\n\n', mode);
fprintf(fid, '| table | level | n | m^3 | H/delta | eta1 | paper kappa | repo kappa | paper it | repo it | relres | status | notes |\n');
fprintf(fid, '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n');
for k = 1:numel(results)
    r = results(k);
    fprintf(fid, '| %d | %s | %d | %d | %.4g | %.3g | %.2f | %.2f | %d | %d | %.2e | %s | %s |\n', ...
        r.table, r.level, r.n, r.m^3, r.Hdelta, r.eta1, r.paperKappa, r.repoKappa, ...
        r.paperIter, r.repoIter, r.relres, r.status, r.notes);
end
fprintf(fid, '\nLarge cells are blocked by default unless `opts.allowLarge=true` is supplied and the memory estimate is below the configured gate.\n');
end
