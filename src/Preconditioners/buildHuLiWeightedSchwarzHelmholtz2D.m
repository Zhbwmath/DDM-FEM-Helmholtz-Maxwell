function method = buildHuLiWeightedSchwarzHelmholtz2D(node, elem, bdFlag, k, parts, opts)
% BUILDHULIWEIGHTEDSCHWARZHELMHOLTZ2D  Hu-Li harmonic coarse space and WASI hybrid.

if nargin < 6 || isempty(opts), opts = struct(); end
totalTimer = tic;
opts = localOptions(opts);
pde = configurePDE(k, opts);
kappaRef = referenceWaveNumber(pde, opts);
fine = setupFineSpace(node, elem, bdFlag, pde, opts.degree);
energyFactorTime = 0;
energySolve = opts.energySolve;
if opts.cacheEnergySolver || opts.cacheEnergyAdjoint
    if isempty(energySolve)
        energyFactorTimer = tic;
        energySolve = energySolverHandle(fine.energy);
        energyFactorTime = toc(energyFactorTimer);
    elseif ~isa(energySolve, 'function_handle')
        error('buildHuLiWeightedSchwarzHelmholtz2D:energySolve', ...
            'opts.energySolve must be a function handle.');
    end
end
if opts.cacheEnergySolver
    fine.energySolve = energySolve;
end
[partData, pouError] = preparePartData(fine, parts);
[solverMode, storedLuGB] = chooseSolverMode(partData, opts);

nSub = numel(partData);
localData = cell(nSub, 1);
if opts.useParfor
    parfor s = 1:nSub
        localData{s} = setupOneLocal(fine, pde, partData(s), ...
            kappaRef, opts, solverMode);
    end
else
    for s = 1:nSub
        localData{s} = setupOneLocal(fine, pde, partData(s), ...
            kappaRef, opts, solverMode);
    end
end

globalBasisTimer = tic;
[Zraw, localInfo, localSolvers] = assembleGlobalBasis(localData, fine.N);
globalBasisTime = toc(globalBasisTimer);
clear localData

rankTimer = tic;
[Z, rankInfo] = independentBasis(Zraw, opts.rankTolerance, opts.rankMethod);
rankTime = toc(rankTimer);

coarseAssemblyTimer = tic;
AH = Z' * fine.A * Z;
coarseAssemblyTime = toc(coarseAssemblyTimer);

coarseFactorTimer = tic;
[coarseSolve, coarseSolveAdjoint] = coarseSolvers(AH, opts.coarseSolverMode);
coarseFactorTime = toc(coarseFactorTimer);

energyAdjointTrial = [];
energyAdjointTime = 0;
if opts.cacheEnergyAdjoint && ~isempty(Z)
    energyAdjointTimer = tic;
    energyAdjointTrial = buildEnergyAdjointTrial( ...
        energySolve, fine.A, Z, opts.energyAdjointBlockSize);
    energyAdjointTime = toc(energyAdjointTimer);
end

    function y = applyM0Inverse(r)
        if isempty(Z)
            y = zeros(fine.N, size(r, 2));
        else
            y = Z * coarseSolve(Z' * r);
        end
    end

    function y = applyWASIInverse(r)
        y = zeros(fine.N, size(r, 2));
        for j = 1:nSub
            data = localSolvers{j};
            if isempty(data.idx), continue; end
            zj = solveFactor(data.solver, r(data.idx, :));
            y(data.idx, :) = y(data.idx, :) + data.weight .* zj;
        end
    end

    function y = applyResidual(r)
        z0 = applyM0Inverse(r);
        zw = applyWASIInverse(r);
        y = z0 + zw - applyM0Inverse(fine.A * zw);
    end

    function y = apply(v)
        y = applyResidual(fine.A * v);
    end

coarseSpace = struct();
coarseSpace.trial = Z;
coarseSpace.test = Z;
coarseSpace.AH = AH;
coarseSpace.solve = coarseSolve;
coarseSpace.solveAdjoint = coarseSolveAdjoint;
coarseSpace.description = sprintf('Hu-Li %s Helmholtz-harmonic coarse space', ...
    lower(opts.coarseType));
coarseSpace.coarseType = lower(opts.coarseType);
coarseSpace.kappaRef = kappaRef;
coarseSpace.rawDimension = size(Zraw, 2);
coarseSpace.dimension = size(Z, 2);
coarseSpace.rankInfo = rankInfo;
coarseSpace.energyAdjointTrial = energyAdjointTrial;

localSizes = arrayfun(@(x) x.localDof, localInfo);
selectedSizes = arrayfun(@(x) x.selectedDimension, localInfo);
method = struct();
method.apply = @apply;
method.applyResidual = @applyResidual;
method.applyM0Inverse = @applyM0Inverse;
method.applyWASIInverse = @applyWASIInverse;
method.A = fine.A;
method.energy = fine.energy;
method.pde = pde;
method.kappaRef = kappaRef;
method.degree = fine.degree;
method.coarseType = lower(opts.coarseType);
method.fineSpace = fine;
method.coarseSpace = coarseSpace;
method.local = localInfo;
method.stats = struct('nSubdomains', nSub, ...
    'localDofMin', min(localSizes), 'localDofMax', max(localSizes), ...
    'localDofMean', mean(localSizes), ...
    'selectedMin', min(selectedSizes), 'selectedMax', max(selectedSizes), ...
    'selectedMean', mean(selectedSizes), ...
    'rawCoarseDimension', size(Zraw, 2), ...
    'coarseDimension', size(Z, 2), ...
    'coarseToMaxLocalRatio', size(Z, 2) / max(localSizes), ...
    'partitionUnityError', pouError, ...
    'solverModeRequested', opts.solverMode, ...
    'solverModeEffective', solverMode, ...
    'estimatedStoredLuGB', storedLuGB);
method.timing = struct('totalSetup', toc(totalTimer), ...
    'localAssemblySum', sum([localInfo.assemblyTime]), ...
    'localFactorSum', sum([localInfo.factorTime]), ...
    'harmonicExtensionSum', sum([localInfo.harmonicTime]), ...
    'coarseSelectionSum', sum([localInfo.coarseSelectionTime]), ...
    'globalBasisAssembly', globalBasisTime, ...
    'rankReduction', rankTime, ...
    'coarseAssembly', coarseAssemblyTime, ...
    'coarseFactorization', coarseFactorTime, ...
    'energyFactorization', energyFactorTime, ...
    'energyAdjointCache', energyAdjointTime);
method.options = opts;
end


function opts = localOptions(opts)
defaults = struct();
defaults.degree = 1;
defaults.coarseType = 'spectral';
defaults.rho = [];
defaults.nu = [];
defaults.beta = [];
defaults.kappaRef = [];
defaults.epsilon = [];
defaults.eta = [];
defaults.solverMode = 'adaptive';
defaults.useParfor = false;
defaults.localStoredLuLimitGB = 200;
defaults.localLuFillConstant = 20;
defaults.denseEigenLimit = 256;
defaults.initialEigenCount = 24;
defaults.eigenTolerance = 1e-10;
defaults.eigenMaxIterations = 1000;
defaults.rankTolerance = 1e-10;
defaults.rankMethod = 'none';
defaults.coarseSolverMode = 'lu';
defaults.energySolve = [];
defaults.cacheEnergySolver = false;
defaults.cacheEnergyAdjoint = false;
defaults.energyAdjointBlockSize = 64;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end

if ~ismember(opts.degree, [1, 2, 3])
    error('buildHuLiWeightedSchwarzHelmholtz2D:degree', ...
        'Only P1, P2, and P3 fine spaces are supported.');
end
if strcmpi(opts.coarseType, 'spectral') && isempty(opts.rho)
    error('buildHuLiWeightedSchwarzHelmholtz2D:rho', ...
        'The spectral coarse space requires opts.rho.');
end
if strcmpi(opts.coarseType, 'economic') && isempty(opts.nu) && isempty(opts.beta)
    error('buildHuLiWeightedSchwarzHelmholtz2D:nu', ...
        'The economic coarse space requires opts.nu or opts.beta.');
end
if ~ismember(lower(opts.rankMethod), {'none', 'gram', 'qr'})
    error('buildHuLiWeightedSchwarzHelmholtz2D:rankMethod', ...
        'rankMethod must be "none", "gram", or "qr".');
end
if ~ismember(lower(opts.coarseSolverMode), {'lu', 'storedlu', 'direct', 'backslash'})
    error('buildHuLiWeightedSchwarzHelmholtz2D:coarseSolverMode', ...
        'coarseSolverMode must be "lu" or "direct".');
end
end


function pde = configurePDE(k, opts)
pde = normalizeHelmholtzPDE(k);
if ~isempty(opts.epsilon)
    pde.epsilon = opts.epsilon;
end
if ~isempty(opts.eta)
    pde.eta = opts.eta;
end
end


function kappaRef = referenceWaveNumber(pde, opts)
if ~isempty(opts.kappaRef)
    kappaRef = opts.kappaRef;
elseif isnumeric(pde.k) && isscalar(pde.k)
    kappaRef = abs(pde.k);
else
    error('buildHuLiWeightedSchwarzHelmholtz2D:kappaRef', ...
        'Variable wave-number input requires opts.kappaRef.');
end
if ~(isnumeric(kappaRef) && isscalar(kappaRef) && isfinite(kappaRef) && kappaRef > 0)
    error('buildHuLiWeightedSchwarzHelmholtz2D:kappaRef', ...
        'opts.kappaRef must be a positive finite scalar.');
end
end


function fine = setupFineSpace(node, elem, bdFlag, pde, degree)
baseElem = elem(:, 1:3);
baseNode = node(1:max(baseElem(:)), :);
if degree == 1
    fineNode = baseNode;
    fineElem = baseElem;
    p1ToFine = speye(size(baseNode, 1));
else
    if size(elem, 2) == 3
        [fineNode, fineElem] = extendMesh2D(baseNode, baseElem, degree);
    else
        fineNode = node;
        fineElem = elem;
    end
    if degree == 2
        p1ToFine = prolongate_P1_P2(baseNode, baseElem);
    else
        p1ToFine = prolongate_P1_P3(baseNode, baseElem);
    end
end

K = assembleStiffness2D(fineNode, fineElem, degree);
M = assembleMass2D(fineNode, fineElem, degree);
Mb = assembleBoundaryMass2D(fineNode, fineElem, bdFlag, degree);
A = assembleHelmholtz2D(fineNode, fineElem, bdFlag, pde, [], [], degree);
qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
energy = K + assembleWeightedMass2D(fineNode, fineElem, degree, qfun);

fine = struct();
fine.dim = 2;
fine.form = 'standard';
fine.degree = degree;
fine.node = fineNode;
fine.elem = fineElem;
fine.bdFlag = bdFlag;
fine.baseNode = baseNode;
fine.baseElem = baseElem;
fine.baseBdFlag = bdFlag;
fine.K = K;
fine.M = M;
fine.boundaryMass = Mb;
fine.A = A;
fine.energy = energy;
fine.pde = pde;
fine.helmholtzInput = pde;
fine.p1ToFine = p1ToFine;
fine.baseToFine = p1ToFine;
fine.N = size(fineNode, 1);
end


function [partData, pouError] = preparePartData(fine, parts)
nSub = numel(parts);
partData = repmat(struct('elemIdx', [], 'idx', [], 'rawWeight', [], ...
    'weight', []), nSub, 1);
weightSum = zeros(fine.N, 1);
covered = false(fine.N, 1);
useWeightFun = isfield(parts, 'weightFun') && ~isempty(parts(1).weightFun);

for s = 1:nSub
    eIdx = parts(s).elemIdx(:);
    idx = unique(fine.elem(eIdx, :));
    if useWeightFun
        raw = max(parts(s).weightFun(fine.node(idx, 1), fine.node(idx, 2)), 0);
    else
        raw = ones(numel(idx), 1);
    end
    partData(s).elemIdx = eIdx;
    partData(s).idx = idx(:);
    partData(s).rawWeight = raw(:);
    weightSum(idx) = weightSum(idx) + raw(:);
    covered(idx) = true;
end

if any(~covered)
    error('buildHuLiWeightedSchwarzHelmholtz2D:partitionCoverage', ...
        'The supplied subdomains do not cover every fine-space DOF.');
end
if any(weightSum <= 0)
    error('buildHuLiWeightedSchwarzHelmholtz2D:partitionWeights', ...
        'Partition weights must have a positive accumulated value at every DOF.');
end

normalizedSum = zeros(fine.N, 1);
for s = 1:nSub
    idx = partData(s).idx;
    partData(s).weight = partData(s).rawWeight ./ weightSum(idx);
    normalizedSum(idx) = normalizedSum(idx) + partData(s).weight;
end
pouError = norm(normalizedSum - 1, inf);
end


function [mode, estimatedGB] = chooseSolverMode(partData, opts)
sizes = arrayfun(@(p) numel(p.idx), partData);
bytes = sum(16 * opts.localLuFillConstant .* sizes .* log2(max(sizes, 2)));
estimatedGB = bytes / 2^30;
if strcmpi(opts.solverMode, 'adaptive')
    if estimatedGB <= opts.localStoredLuLimitGB
        mode = 'lu';
    else
        mode = 'direct';
    end
else
    mode = lower(opts.solverMode);
end
if ~ismember(mode, {'lu', 'direct', 'backslash'})
    error('buildHuLiWeightedSchwarzHelmholtz2D:solverMode', ...
        'solverMode must be "adaptive", "lu", or "direct".');
end
if strcmp(mode, 'backslash'), mode = 'direct'; end
end


function data = setupOneLocal(fine, pde, part, kappaRef, opts, solverMode)
assemblyTimer = tic;
idx = part.idx;
eIdx = part.elemIdx;
g2l = zeros(fine.N, 1);
g2l(idx) = (1:numel(idx))';
localNode = fine.node(idx, :);
localElem = g2l(fine.elem(eIdx, :));
localBdFlag = localBoundaryFlags(localElem(:, 1:3));

A = assembleHelmholtz2D(localNode, localElem, localBdFlag, ...
    pde, [], [], fine.degree);
K = assembleStiffness2D(localNode, localElem, fine.degree);
qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
D = K + assembleWeightedMass2D(localNode, localElem, fine.degree, qfun);
Mb = assembleBoundaryMass2D(localNode, localElem, localBdFlag, fine.degree);
boundaryIdx = find(sum(abs(Mb), 2) > 0);
rhs = Mb(:, boundaryIdx);
assemblyTime = toc(assemblyTimer);

factorTimer = tic;
solver = factorMatrix(A, solverMode);
factorTime = toc(factorTimer);
coarseTimer = tic;
switch lower(opts.coarseType)
    case 'spectral'
        harmonicTimer = tic;
        H = solveFactor(solver, rhs);
        harmonicTime = toc(harmonicTimer);
        harmonicResidual = norm(A * H - rhs, 'fro') / max(1, norm(rhs, 'fro'));
        [modes, spectral] = spectralLocalModes(H, D, part.weight, opts);
        harmonicDimension = size(H, 2);
        traceDimension = numel(boundaryIdx);
        tracePeriodicityError = NaN;
    case 'economic'
        nu = economicResolution(kappaRef, opts);
        [traceBasis, tracePeriodicityError] = periodicQuadraticTrace( ...
            localNode(boundaryIdx, :), nu);
        economicRhs = rhs * traceBasis;
        harmonicTimer = tic;
        modes = solveFactor(solver, economicRhs);
        harmonicTime = toc(harmonicTimer);
        harmonicResidual = norm(A * modes - economicRhs, 'fro') / ...
            max(1, norm(economicRhs, 'fro'));
        harmonicDimension = size(modes, 2);
        traceDimension = size(traceBasis, 2);
        spectral = emptySpectralInfo();
        spectral.selectedDimension = size(modes, 2);
    otherwise
        error('buildHuLiWeightedSchwarzHelmholtz2D:coarseType', ...
            'Unknown coarseType "%s".', opts.coarseType);
end
coarseSelectionTime = toc(coarseTimer);

weightedModes = part.weight .* modes;
data = struct();
data.idx = idx;
data.weight = part.weight;
data.solver = solver;
data.weightedModes = weightedModes;
data.info = struct('localDof', numel(idx), ...
    'boundaryDof', numel(boundaryIdx), ...
    'harmonicDimension', harmonicDimension, ...
    'selectedDimension', size(modes, 2), ...
    'traceDimension', traceDimension, ...
    'harmonicResidual', harmonicResidual, ...
    'tracePeriodicityError', tracePeriodicityError, ...
    'eigenvalues', spectral.eigenvalues, ...
    'selectedEigenvalues', spectral.selectedEigenvalues, ...
    'firstExcludedEigenvalue', spectral.firstExcludedEigenvalue, ...
    'energyOrthonormalityError', spectral.energyOrthonormalityError, ...
    'eigenResidual', spectral.eigenResidual, ...
    'hermitianLeftError', spectral.hermitianLeftError, ...
    'hermitianRightError', spectral.hermitianRightError, ...
    'assemblyTime', assemblyTime, 'factorTime', factorTime, ...
    'harmonicTime', harmonicTime, ...
    'coarseSelectionTime', coarseSelectionTime);
end


function info = emptySpectralInfo()
info = struct('eigenvalues', [], 'selectedEigenvalues', [], ...
    'firstExcludedEigenvalue', NaN, 'selectedDimension', 0, ...
    'energyOrthonormalityError', NaN, 'eigenResidual', NaN, ...
    'hermitianLeftError', NaN, 'hermitianRightError', NaN);
end


function [modes, info] = spectralLocalModes(H, D, weight, opts)
XH = weight .* H;
S = H' * (weight .* (D * XH));
T = H' * (D * H);
leftError = norm(S - S', 'fro') / max(1, norm(S, 'fro'));
rightError = norm(T - T', 'fro') / max(1, norm(T, 'fro'));
S = (S + S') / 2;
T = (T + T') / 2;
info = emptySpectralInfo();
info.hermitianLeftError = leftError;
info.hermitianRightError = rightError;

[V, lambda] = generalizedHermitianEigenpairs(S, T, opts, opts.rho^2);
lambda = real(lambda(:));
[lambda, order] = sort(lambda, 'descend');
V = V(:, order);
selected = lambda >= opts.rho^2 * (1 - 10 * opts.eigenTolerance);
Vsel = V(:, selected);
modes = H * Vsel;

info.eigenvalues = lambda;
info.selectedEigenvalues = lambda(selected);
info.selectedDimension = nnz(selected);
excluded = lambda(~selected);
if ~isempty(excluded)
    info.firstExcludedEigenvalue = max(excluded);
end
if isempty(Vsel)
    info.energyOrthonormalityError = 0;
    info.eigenResidual = 0;
else
    info.energyOrthonormalityError = norm(Vsel' * T * Vsel - eye(size(Vsel, 2)), 'fro');
    residual = S * Vsel - T * (Vsel .* reshape(lambda(selected).', 1, []));
    info.eigenResidual = norm(residual, 'fro') / max(1, norm(S * Vsel, 'fro'));
end
end


function [V, lambda] = generalizedHermitianEigenpairs(S, T, opts, threshold)
n = size(S, 1);
if n == 0
    V = zeros(0, 0);
    lambda = zeros(0, 1);
    return;
end

if n == 1 || n <= opts.denseEigenLimit
    [V, lambda] = eig(full(S), full(T), 'vector');
    good = isfinite(lambda) & abs(imag(lambda)) <= 100 * opts.eigenTolerance;
    V = V(:, good);
    lambda = real(lambda(good));
    return;
end

S = sparse((S + S') / 2);
T = positiveDefiniteSparseMass((T + T') / 2, opts);
nev = min(max(1, opts.initialEigenCount), n - 1);
eigsOpts = struct('tol', opts.eigenTolerance, ...
    'maxit', opts.eigenMaxIterations);
while true
    try
        [V, L] = eigs(S, T, nev, 'largestreal', eigsOpts);
    catch
        [V, L] = eigs(S, T, nev, 'lm', eigsOpts);
    end
    lambda = real(diag(L));
    [lambda, order] = sort(lambda, 'descend');
    V = V(:, order);
    if min(lambda) < threshold || nev == n - 1
        break;
    end
    nev = min(n - 1, max(nev + 1, 2 * nev));
end

if nev == n - 1 && min(lambda) >= threshold
    try
        [vs, ls] = eigs(S, T, 1, 'smallestreal', eigsOpts);
    catch
        [vs, ls] = eigs(S, T, 1, 'sm', eigsOpts);
    end
    V = [V, vs];
    lambda = [lambda; real(ls(1,1))];
end
end


function T = positiveDefiniteSparseMass(T, opts)
T = sparse(T);
[~, flag] = chol(T);
if flag == 0
    return;
end
scale = max(1, norm(T, 1));
shift = max(opts.eigenTolerance, eps) * scale;
for attempt = 1:8
    Ttry = T + shift * speye(size(T, 1));
    [~, flag] = chol(Ttry);
    if flag == 0
        T = Ttry;
        return;
    end
    shift = 10 * shift;
end
error('buildHuLiWeightedSchwarzHelmholtz2D:eigenMassMatrix', ...
    'Spectral generalized eigenproblem mass matrix is not positive definite.');
end


function nu = economicResolution(kappaRef, opts)
if ~isempty(opts.nu)
    nu = round(opts.nu);
else
    nu = round(kappaRef^(1 - opts.beta));
end
nu = max(1, nu);
end


function [B, periodicityError] = periodicQuadraticTrace(boundaryNode, nu)
t = rectangularBoundaryParameter(boundaryNode);
B = evalPeriodicQuadraticBasis(t, nu);
B0 = evalPeriodicQuadraticBasis(0, nu);
B1 = evalPeriodicQuadraticBasis(1, nu);
periodicityError = norm(B0 - B1, inf);
end


function t = rectangularBoundaryParameter(x)
xmin = min(x(:,1)); xmax = max(x(:,1));
ymin = min(x(:,2)); ymax = max(x(:,2));
Lx = xmax - xmin; Ly = ymax - ymin;
perimeter = 2 * (Lx + Ly);
tol = max(1, perimeter) * 1e-10;
if Lx <= tol || Ly <= tol
    error('buildHuLiWeightedSchwarzHelmholtz2D:economicBoundary', ...
        'Economic traces require a nondegenerate rectangular subdomain.');
end

s = nan(size(x,1), 1);
bottom = abs(x(:,2) - ymin) <= tol;
s(bottom) = x(bottom,1) - xmin;
right = isnan(s) & abs(x(:,1) - xmax) <= tol;
s(right) = Lx + x(right,2) - ymin;
top = isnan(s) & abs(x(:,2) - ymax) <= tol;
s(top) = Lx + Ly + xmax - x(top,1);
left = isnan(s) & abs(x(:,1) - xmin) <= tol;
s(left) = 2 * Lx + Ly + ymax - x(left,2);
if any(isnan(s))
    error('buildHuLiWeightedSchwarzHelmholtz2D:economicBoundary', ...
        'Economic traces currently require axis-aligned rectangular subdomains.');
end
t = mod(s / perimeter, 1);
end


function B = evalPeriodicQuadraticBasis(t, nu)
t = mod(t(:), 1);
nDof = 2 * nu;
u = t * nu;
e = floor(u);
xi = u - e;
phi = lagrange1D(2, xi);
left = 2 * e + 1;
mid = 2 * e + 2;
right = mod(2 * e + 2, nDof) + 1;
rows = repmat((1:numel(t))', 1, 3);
cols = [left, mid, right];
B = sparse(rows(:), cols(:), phi(:), numel(t), nDof);
end


function [Zraw, info, solvers] = assembleGlobalBasis(localData, N)
nSub = numel(localData);
blocks = cell(nSub, 1);
info = repmat(localData{1}.info, nSub, 1);
solvers = cell(nSub, 1);
for s = 1:nSub
    data = localData{s};
    nMode = size(data.weightedModes, 2);
    [ii, jj] = ndgrid(data.idx, 1:nMode);
    blocks{s} = sparse(ii(:), jj(:), data.weightedModes(:), N, nMode);
    info(s) = data.info;
    solvers{s} = struct('idx', data.idx, 'weight', data.weight, ...
        'solver', data.solver);
end
Zraw = [blocks{:}];
end


function [Z, info] = independentBasis(Zraw, tolerance, method)
rawDimension = size(Zraw, 2);
if rawDimension == 0
    Z = sparse(size(Zraw, 1), 0);
    info = struct('rawDimension', 0, 'dimension', 0, ...
        'relativeTolerance', tolerance, 'method', lower(method), ...
        'pivot', zeros(1,0), ...
        'diagonalR', zeros(0,1));
    return;
end

columnNorm = sqrt(full(sum(abs(Zraw).^2, 1)));
nonzeroColumn = columnNorm > tolerance * max(max(columnNorm), eps);
active = find(nonzeroColumn);
if isempty(active)
    Z = sparse(size(Zraw, 1), 0);
    info = struct('rawDimension', rawDimension, 'dimension', 0, ...
        'relativeTolerance', tolerance, 'method', lower(method), ...
        'pivot', zeros(1,0), 'diagonalR', zeros(0,1));
    return;
end

scaled = Zraw(:, active) * spdiags(1 ./ columnNorm(active).', ...
    0, numel(active), numel(active));
switch lower(method)
    case 'none'
        Z = scaled;
        rankZ = size(scaled, 2);
        localPivot = 1:rankZ;
        d = columnNorm(active).';
        selectedLocal = localPivot;
    case 'gram'
        gram = full(scaled' * scaled);
        gram = (gram + gram') / 2;
        [~, R, localPivot] = qr(gram, 'vector');
        d = abs(diag(R));
        scale = max(d);
        rankZ = nnz(d > tolerance^2 * max(scale, eps));
        selectedLocal = localPivot(1:rankZ);
        Z = scaled(:, selectedLocal);
    case 'qr'
        [Q, R, localPivot] = qr(scaled, 0);
        d = abs(diag(R));
        scale = max(d);
        rankZ = nnz(d > tolerance * max(scale, eps));
        selectedLocal = localPivot(1:rankZ);
        Z = Q(:, 1:rankZ);
end
pivot = active(localPivot);
info = struct('rawDimension', rawDimension, 'dimension', rankZ, ...
    'relativeTolerance', tolerance, 'method', lower(method), ...
    'pivot', pivot, 'selectedColumns', active(selectedLocal), ...
    'diagonalR', d);
end


function [solve, solveAdjoint] = coarseSolvers(A, mode)
if isempty(A)
    solve = @(r) zeros(0, size(r,2));
    solveAdjoint = solve;
    return;
end
switch lower(mode)
    case {'direct', 'backslash'}
        solve = @(r) A \ r;
        solveAdjoint = @(r) A' \ r;
    otherwise
        factor = factorMatrix(A, 'lu');
        solve = @(r) solveFactor(factor, r);
        solveAdjoint = @(r) solveFactorAdjoint(factor, r);
end
end


function trial = buildEnergyAdjointTrial(solveEnergy, A, Z, blockSize)
nCoarse = size(Z, 2);
trial = complex(zeros(size(Z, 1), nCoarse));
for first = 1:blockSize:nCoarse
    columns = first:min(first + blockSize - 1, nCoarse);
    trial(:, columns) = solveEnergy(A' * Z(:, columns));
end
end


function solve = energySolverHandle(D)
try
    R = chol(D);
    solve = @(b) R \ (R' \ b);
catch
    factor = factorMatrix(D, 'lu');
    solve = @(b) solveFactor(factor, b);
end
end


function solver = factorMatrix(A, mode)
switch lower(mode)
    case {'direct', 'backslash'}
        solver = struct('mode', 'direct', 'A', A);
    case {'lu', 'storedlu'}
        [L, U, p, q] = lu(A, 'vector');
        solver = struct('mode', 'lu', 'L', L, 'U', U, ...
            'p', p(:), 'q', q(:));
    otherwise
        error('buildHuLiWeightedSchwarzHelmholtz2D:factorMode', ...
            'Unknown factorization mode "%s".', mode);
end
end


function x = solveFactor(solver, b)
if strcmp(solver.mode, 'direct')
    x = solver.A \ b;
else
    x = zeros(size(b));
    x(solver.q, :) = solver.U \ (solver.L \ b(solver.p, :));
end
end


function x = solveFactorAdjoint(solver, b)
if strcmp(solver.mode, 'direct')
    x = solver.A' \ b;
else
    x = zeros(size(b));
    x(solver.p, :) = solver.L' \ (solver.U' \ b(solver.q, :));
end
end


function bdFlag = localBoundaryFlags(vertexElem)
edgePairs = [2 3; 3 1; 1 2];
allEdges = [vertexElem(:, edgePairs(1,:)); ...
    vertexElem(:, edgePairs(2,:)); vertexElem(:, edgePairs(3,:))];
[~, ~, edgeId] = unique(sort(allEdges, 2), 'rows');
counts = accumarray(edgeId, 1);
bdFlag = reshape(counts(edgeId) == 1, size(vertexElem, 1), 3);
end
