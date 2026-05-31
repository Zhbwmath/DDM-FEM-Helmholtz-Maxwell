% VERIFY_GSZ20_ORAS_IMP_HELMHOLTZ  GSZ20 impedance ORAS table reproduction.

fprintf('========== GSZ20 impedance ORAS table reproduction ==========\n\n');

cfg = defaultConfig();
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
outDir = fullfile(repoRoot, 'docs', 'GSZ20_oras_imp_helmholtz', 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

results = [];
writeResultFiles(results, outDir, cfg);

fprintf('Estimated largest Table 1-3 case memory: %.1f GB\n\n', estimateLargestMemoryGB(cfg));

fprintf('Table 1(a): eps_prob=eps_prec=k^(1+alpha+0.1)\n');
for k = cfg.tableK
    for alpha = cfg.alphaValues
        beta = alpha + 0.1;
        epsv = k^(1 + beta);
        results = runAndCheckpoint(results, outDir, 'table1a', k, alpha, [], epsv, epsv, ...
            'impedance', 'random', cfg);
    end
end

fprintf('\nTable 1(b): eps_prob=eps_prec=k\n');
for k = cfg.tableK
    for alpha = cfg.alphaValues
        results = runAndCheckpoint(results, outDir, 'table1b', k, alpha, [], k, k, ...
            'impedance', 'random', cfg);
    end
end

fprintf('\nTable 2: eps_prob=0, eps_prec in {k,0}\n');
for epsMode = ["k", "0"]
    for k = cfg.tableK
        for alpha = cfg.alphaValues
            epsPrec = k * strcmp(epsMode, "k");
            results = runAndCheckpoint(results, outDir, "table2_epsprec_" + epsMode, ...
                k, alpha, [], 0, epsPrec, 'impedance', 'random', cfg);
        end
    end
end

fprintf('\nTable 3: fixed M, eps_prob=0\n');
for panel = ["random_epsprec_k", "random_epsprec_0", "zero_epsprec_k"]
    for k = cfg.tableK
        for mFixed = cfg.fixedMValues
            epsPrec = k;
            x0Type = 'random';
            if panel == "random_epsprec_0"
                epsPrec = 0;
            elseif panel == "zero_epsprec_k"
                x0Type = 'zero';
            end
            results = runAndCheckpoint(results, outDir, "table3_" + panel, k, [], ...
                mFixed, 0, epsPrec, 'impedance', x0Type, cfg);
        end
    end
end

fprintf('\nTable 4: skipped because strict Dirichlet setup is ambiguous for this RHS.\n');

writeReport(results, fullfile(repoRoot, 'docs', 'GSZ20_oras_imp_helmholtz', ...
    'graham_spence_zou_tables_1_4_codex_note.md'), cfg);

fprintf('\nWrote GSZ20 table results to %s\n', outDir);
fprintf('========== GSZ20 table reproduction complete ==========\n');


function cfg = defaultConfig()
cfg.seed = 20240528;
cfg.tol = 1e-6;
cfg.maxIter = 100;
cfg.gmresBasis = 100;
cfg.sparseLUFillConstant = 40;
cfg.rounding = 'round';
cfg.etaConvention = 'eta=k';
cfg.tableK = [40, 60, 80, 100, 120, 140];
cfg.alphaValues = [0.2, 0.3, 0.4, 0.5];
cfg.fixedMValues = [4, 8, 16];
cfg.nfineCap = inf;
capText = getenv('GSZ20_NFINE_CAP');
if ~isempty(capText)
    cfg.nfineCap = str2double(capText);
end
cfg.date = char(datetime('today', 'Format', 'yyyy-MM-dd'));
end


function gb = estimateLargestMemoryGB(cfg)
k = 140;
mVals = [round(k^0.2), round(k^0.3), round(k^0.4), round(k^0.5), 4, 8, 16];
gb = 0;
for mCoarse = unique(mVals)
    gb = max(gb, estimateCaseMemoryGB(k, mCoarse, cfg));
end
end


function gb = estimateCaseMemoryGB(k, mCoarse, cfg)
nFine = mCoarse * max(1, round(round(k^1.5) / mCoarse));
N = (nFine + 1)^2;
NT = 2 * nFine^2;
step = nFine / mCoarse;
subBytes = 0;
for j = 0:mCoarse
    for i = 0:mCoarse
        nx = min(nFine, (i+1)*step) - max(0, (i-1)*step) + 1;
        ny = min(nFine, (j+1)*step) - max(0, (j-1)*step) + 1;
        nloc = nx * ny;
        subBytes = subBytes + 336*nloc + ...
            16*cfg.sparseLUFillConstant*nloc*log2(max(nloc, 2));
    end
end
gmresBytes = 16*N*(cfg.gmresBasis + 3);
bytes = 112*N + 48*NT + subBytes + gmresBytes;
gb = bytes / 1024^3;
end


function r = runCase(panel, k, alpha, mFixed, epsProb, epsPrec, localBc, x0Type, cfg)
if isempty(mFixed)
    mCoarse = chooseM(k, alpha, cfg.rounding);
else
    mCoarse = mFixed;
end
[nFine, scaled] = chooseNfine(k, mCoarse, cfg.nfineCap);

[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/nFine);
[A, rhs] = assembleProblem(node, elem, bdFlag, k, epsProb);
sub = buildVertexPatchSubdomains(node, elem, nFine, mCoarse, k, epsPrec, localBc);
applyB = @(v) applyPreconditioner(v, sub, size(A, 1));
btilde = applyB(rhs);
C = @(v) applyB(A * v);

n = size(A, 1);
rng(cfg.seed, 'twister');
if strcmpi(x0Type, 'zero')
    x0 = zeros(n, 1);
else
    x0 = rand(n, 1);
end

res0 = norm(btilde - C(x0));
tic;
restart = min([n, cfg.gmresBasis, cfg.maxIter]);
maxCycles = max(1, ceil(cfg.maxIter / restart));
[~, flag, relres, iter, resvec] = gmres(C, btilde, restart, cfg.tol, maxCycles, [], [], x0);
elapsed = toc;
if ~isempty(resvec)
    its = numel(resvec) - 1;
elseif numel(iter) == 2
    its = (iter(1)-1) * restart + iter(2);
else
    its = iter;
end
precRelres = norm(btilde - C(x0)) / max(res0, eps);
if ~isempty(resvec)
    precRelres = resvec(end) / max(resvec(1), eps);
end

r = struct('panel', char(panel), 'k', k, 'alpha', scalarOrNaN(alpha), ...
    'mFixed', scalarOrNaN(mFixed), 'mCoarse', mCoarse, 'nFine', nFine, ...
    'h', 1/nFine, 'nDofs', n, 'epsProb', epsProb, 'epsPrec', epsPrec, ...
    'localBc', localBc, 'x0Type', x0Type, 'seed', cfg.seed, ...
    'gmresIts', its, 'flag', flag, 'relres', relres, ...
    'precRelres', precRelres, 'elapsedSeconds', elapsed, ...
    'scaledNfine', scaled);

fprintf('%-34s k=%-3g alpha=%4.1f M=%-2d N=%-4d dof=%-7d bc=%-9s epsP=%-8.3g its=%-4d flag=%d rel=%.2e%s\n', ...
    r.panel, r.k, r.alpha, r.mCoarse, r.nFine, r.nDofs, r.localBc, ...
    r.epsPrec, r.gmresIts, r.flag, r.precRelres, ternary(r.scaledNfine, ' scaled', ''));
end


function results = runAndCheckpoint(results, outDir, panel, k, alpha, mFixed, epsProb, epsPrec, localBc, x0Type, cfg)
r = runCase(panel, k, alpha, mFixed, epsProb, epsPrec, localBc, x0Type, cfg);
results = [results; r];
writeResultFiles(results, outDir, cfg);
end


function [A, rhs] = assembleProblem(node, elem, bdFlag, k, epsProb)
K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
Mb = assembleBoundaryMass2D(node, elem, bdFlag, 1);
A = K - (k^2 + 1i * epsProb) * M - 1i * k * Mb;

d = [1, 1] / sqrt(2);
uex = exp(1i * k * (node(:,1) * d(1) + node(:,2) * d(2)));
rhs = M * (-1i * epsProb * uex);

bdLoad = boundaryPlaneWaveLoad(node, elem, bdFlag, k, d);
rhs = rhs + bdLoad;
end


function b = boundaryPlaneWaveLoad(node, elem, bdFlag, k, d)
N = size(node, 1);
b = zeros(N, 1);
edgeVertex = [2 3; 3 1; 1 2];
for eLocal = 1:3
    rows = find(bdFlag(:, eLocal) == 1);
    for rr = rows(:)'
        va = elem(rr, edgeVertex(eLocal, 1));
        vb = elem(rr, edgeVertex(eLocal, 2));
        pa = node(va, :);
        pb = node(vb, :);
        edge = pb - pa;
        L = norm(edge);
        mid = 0.5 * (pa + pb);
        normal = outwardNormal(mid);
        ua = exp(1i * k * dot(pa, d));
        ub = exp(1i * k * dot(pb, d));
        ga = 1i * (k * dot(d, normal) - k) * ua;
        gb = 1i * (k * dot(d, normal) - k) * ub;
        b(va) = b(va) + L * (2 * ga + gb) / 6;
        b(vb) = b(vb) + L * (ga + 2 * gb) / 6;
    end
end
end


function n = outwardNormal(p)
tol = 1e-12;
if abs(p(1)) < tol
    n = [-1, 0];
elseif abs(p(1) - 1) < tol
    n = [1, 0];
elseif abs(p(2)) < tol
    n = [0, -1];
elseif abs(p(2) - 1) < tol
    n = [0, 1];
else
    error('Boundary midpoint is not on the unit-square boundary.');
end
end


function sub = buildVertexPatchSubdomains(node, elem, nFine, mCoarse, k, epsPrec, localBc)
step = nFine / mCoarse;
assert(abs(step - round(step)) < 1e-12, 'Nfine must be divisible by Mcoarse.');
step = round(step);
nNodes1D = nFine + 1;
nodeGrid = reshape(1:size(node, 1), nNodes1D, nNodes1D);
sub((mCoarse + 1)^2) = struct('nodes', [], 'weights', [], 'solve', []);
s = 0;
for j = 0:mCoarse
    for i = 0:mCoarse
        s = s + 1;
        ix = max(0, i-1) * step : min(nFine, (i+1) * step);
        iy = max(0, j-1) * step : min(nFine, (j+1) * step);
        gnodes = nodeGrid(ix + 1, iy + 1);
        gnodes = gnodes(:);
        inPatch = false(size(node, 1), 1);
        inPatch(gnodes) = true;
        elemMask = all(inPatch(elem), 2);
        g2l = zeros(size(node, 1), 1);
        g2l(gnodes) = 1:numel(gnodes);
        lelem = g2l(elem(elemMask, :));
        lnode = node(gnodes, :);
        bdLoc = localBoundaryFlags(lelem);
        K = assembleStiffness2D(lnode, lelem, 1);
        M = assembleMass2D(lnode, lelem, 1);
        Mb = assembleBoundaryMass2D(lnode, lelem, bdLoc, 1);
        A = K - (k^2 + 1i * epsPrec) * M;
        if strcmpi(localBc, 'impedance')
            A = A - 1i * k * Mb;
        elseif strcmpi(localBc, 'dirichlet')
            isBd = localBoundaryNodes(lelem);
            solverNodes = find(~isBd);
            A = A(solverNodes, solverNodes);
            gnodes = gnodes(solverNodes);
        else
            error('Unknown local BC: %s', localBc);
        end
        chi = coarseHat(node(gnodes, :), i, j, mCoarse);
        [L, U, p, q] = lu(A, 'vector');
        sub(s).nodes = gnodes;
        sub(s).weights = chi(:);
        sub(s).solve = @(r) localSolve(L, U, p, q, r);
    end
end
end


function bdFlag = localBoundaryFlags(elem)
edgeVertex = [2 3; 3 1; 1 2];
edges = [sort(elem(:, edgeVertex(1,:)), 2); ...
         sort(elem(:, edgeVertex(2,:)), 2); ...
         sort(elem(:, edgeVertex(3,:)), 2)];
[~, ~, eid] = unique(edges, 'rows');
counts = accumarray(eid, 1);
isBdAll = counts(eid) == 1;
nElem = size(elem, 1);
bdFlag = [isBdAll(1:nElem), isBdAll(nElem+1:2*nElem), isBdAll(2*nElem+1:end)];
end


function isBd = localBoundaryNodes(elem)
bdFlag = localBoundaryFlags(elem);
edgeVertex = [2 3; 3 1; 1 2];
isBd = false(max(elem(:)), 1);
for eLocal = 1:3
    rows = bdFlag(:, eLocal) == 1;
    vv = elem(rows, edgeVertex(eLocal, :));
    isBd(vv(:)) = true;
end
end


function z = applyPreconditioner(r, sub, n)
z = zeros(n, 1);
for s = 1:numel(sub)
    idx = sub(s).nodes;
    chi = sub(s).weights;
    zloc = sub(s).solve(chi .* r(idx));
    z(idx) = z(idx) + chi .* zloc;
end
end


function z = localSolve(L, U, p, q, r)
z = zeros(size(r));
z(q) = U \ (L \ r(p));
end


function chi = coarseHat(xy, i, j, mCoarse)
chi = max(1 - abs(mCoarse * xy(:,1) - i), 0) .* ...
      max(1 - abs(mCoarse * xy(:,2) - j), 0);
end


function m = chooseM(k, alpha, mode)
switch mode
    case 'round'
        m = max(1, round(k^alpha));
    case 'ceil'
        m = max(1, ceil(k^alpha));
    case 'floor'
        m = max(1, floor(k^alpha));
    otherwise
        error('Unknown rounding mode: %s', mode);
end
end


function [nFine, scaled] = chooseNfine(k, mCoarse, cap)
nTarget = round(k^1.5);
nFine = mCoarse * max(1, round(nTarget / mCoarse));
scaled = false;
if nFine > cap
    nFine = mCoarse * max(1, floor(cap / mCoarse));
    scaled = true;
end
end


function writeResultFiles(results, outDir, cfg)
if isempty(results)
    T = emptyResultTable();
else
    T = struct2table(results);
end
writetable(T, fullfile(outDir, 'iteration_counts.csv'));
fid = fopen(fullfile(outDir, 'iteration_counts.json'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(results, PrettyPrint=true));
clear cleanup;

writePanel(T, outDir, 'table1', 'table1.csv');
writePanel(T, outDir, 'table2', 'table2.csv');
writePanel(T, outDir, 'table3', 'table3.csv');
table4 = table("skipped", ...
    "Strict local Dirichlet branch is ambiguous with the plane-wave eps_prob=0 boundary-only RHS; rerun only after RHS or local-boundary convention is clarified.", ...
    'VariableNames', {'status', 'reason'});
writetable(table4, fullfile(outDir, 'table4_status.csv'));

meta = table(string(cfg.etaConvention), string(cfg.rounding), cfg.seed, ...
    cfg.tol, cfg.maxIter, cfg.gmresBasis, cfg.sparseLUFillConstant, cfg.nfineCap, ...
    'VariableNames', {'etaConvention', 'rounding', 'seed', 'tol', 'maxIter', 'gmresBasis', 'sparseLUFillConstant', 'nfineCap'});
writetable(meta, fullfile(outDir, 'run_metadata.csv'));
end


function T = emptyResultTable()
T = table(string.empty(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    string.empty(0,1), string.empty(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), ...
    'VariableNames', {'panel','k','alpha','mFixed','mCoarse','nFine','h', ...
    'nDofs','epsProb','epsPrec','localBc','x0Type','seed','gmresIts', ...
    'flag','relres','precRelres','elapsedSeconds','scaledNfine'});
end


function writePanel(T, outDir, prefix, filename)
if height(T) == 0
    writetable(T, fullfile(outDir, filename));
    return;
end
mask = startsWith(string(T.panel), prefix);
writetable(T(mask, :), fullfile(outDir, filename));
end


function writeReport(results, reportPath, cfg)
old = fileread(reportPath);
marker = '## 10. Codex full Table 1--3 reproduction run';
oldMarker = '## 10. Codex reduced experiment run';
pos = [strfind(old, marker), strfind(old, oldMarker)];
if ~isempty(pos)
    old = strtrim(old(1:min(pos)-1));
end
old = regexprep(old, 'Updated: \d{4}-\d{2}-\d{2}', ['Updated: ' cfg.date], 'once');
fid = fopen(reportPath, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n\n', old);
fprintf(fid, '%s\n\n', marker);
fprintf(fid, 'Updated: %s\n\n', cfg.date);
fprintf(fid, 'Verification entry point: `verify/verify_gsz20_oras_imp_helmholtz.m`\n\n');
fprintf(fid, 'Main utilities: `squaremesh`, `assembleStiffness2D`, `assembleMass2D`, `assembleBoundaryMass2D`, MATLAB `gmres`.\n\n');
fprintf(fid, 'This run intentionally skips internet search and uses the formulation recorded above. ');
fprintf(fid, 'It reproduces the full Table 1--3 parameter grids requested in the note. ');
fprintf(fid, 'Table 4 is skipped because the strict Dirichlet interpretation made the left-preconditioned system degenerate for the boundary-only plane-wave RHS in the earlier diagnostic run.\n\n');
if isinf(cfg.nfineCap)
    fprintf(fid, '- Fine-mesh convention: `Ntarget = round(k^(3/2))`, adjusted to be divisible by `Mcoarse`; no `Nfine` cap is active.\n');
else
    fprintf(fid, '- Fine-mesh convention: `Ntarget = round(k^(3/2))`, adjusted to be divisible by `Mcoarse`; debug cap `Nfine <= %g` is active.\n', cfg.nfineCap);
end
fprintf(fid, '- Coarse convention: `Mcoarse = round(k^alpha)` unless Table 3 fixes `M`.\n');
fprintf(fid, '- Partition: vertex-patch subdomains with tensor-product coarse hat weights; both restriction and prolongation are weighted by `chi_l` as specified in the note.\n');
fprintf(fid, '- Global operator: `S - (k^2 + i eps_prob) M - i k N_boundary`, with plane-wave boundary data.\n');
fprintf(fid, '- Local impedance operator: `S_l - (k^2 + i eps_prec) M_l - i k N_l_boundary`.\n');
fprintf(fid, '- GMRES: left-preconditioned MATLAB `gmres`, tolerance `%.1e`, seed `%d`, `eta=k`, restart basis `%d`, maximum iterations `%d`.\n\n', cfg.tol, cfg.seed, cfg.gmresBasis, cfg.maxIter);
fprintf(fid, 'Estimated largest Table 1--3 memory from the project sparse-LU rule of thumb is %.1f GB with `c_lu=%g`, below the 200 GB permission gate. Actual sparse LU fill can vary with ordering, so large cases should still be monitored.\n\n', estimateLargestMemoryGB(cfg), cfg.sparseLUFillConstant);

T = struct2table(results);
nFlag = sum(T.flag ~= 0);
fprintf(fid, 'Run status: `%d` Table 1--3 cases completed; `%d` converged with `flag=0`; `%d` hit the 100-iteration GMRES cap. The capped rows are retained in the tables with `flag=1` rather than rerun with a larger maximum because this run follows the requested `maxIter=100` setting.\n\n', ...
    height(T), height(T) - nFlag, nFlag);
fprintf(fid, '| panel | k | alpha | M | Nfine | dofs | eps_prob | eps_prec | local_bc | x0 | its | flag | relres | notes |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---|\n');
for i = 1:height(T)
    notes = '';
    if T.scaledNfine(i)
        notes = 'scaled Nfine cap';
    end
    fprintf(fid, '| %s | %g | %g | %d | %d | %d | %.6g | %.6g | %s | %s | %d | %d | %.3g | %s |\n', ...
        T.panel{i}, T.k(i), T.alpha(i), T.mCoarse(i), T.nFine(i), T.nDofs(i), ...
        T.epsProb(i), T.epsPrec(i), T.localBc{i}, T.x0Type{i}, ...
        T.gmresIts(i), T.flag(i), T.precRelres(i), notes);
end
fprintf(fid, '\nResult files are under `docs/GSZ20_oras_imp_helmholtz/results/`.\n');
clear cleanup;
end


function y = scalarOrNaN(x)
if isempty(x)
    y = NaN;
else
    y = x;
end
end


function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
