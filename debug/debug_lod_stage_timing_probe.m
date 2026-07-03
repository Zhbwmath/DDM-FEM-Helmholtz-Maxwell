% DEBUG_LOD_STAGE_TIMING_PROBE  Time internal Helmholtz LOD build stages.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

k = envNumber('LOD_STAGE_PROBE_K', 32);
ell = envNumber('LOD_STAGE_PROBE_ELL', 2);
workers = envNumber('LOD_STAGE_PROBE_WORKERS', 4);
useParfor = logical(envNumber('LOD_STAGE_PROBE_PARFOR', 1));
Ch = envNumber('LOD_STAGE_PROBE_CH', 1);

degree = 1;
HInv = max(1, round(k));
rawFine = ceil(Ch * k^((2 * degree + 1) / (2 * degree)));
hInv = ceil(rawFine / HInv) * HInv;
h = 1 / hInv;
H = 1 / HInv;

fprintf('========== LOD internal stage timing probe ==========\n');
fprintf('k=%g, h^{-1}=%d, H^{-1}=%d, ell=%d, useParfor=%d, workers=%d\n', ...
    k, hInv, HInv, ell, useParfor, workers);
printMemory('start');

if useParfor
    pool = ensureExactParpool(workers);
    fprintf('Active parpool workers confirmed: p.NumWorkers=%d\n', pool.NumWorkers);
end

t = tic;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], h);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], H);
stage('mesh', t);
fprintf('N=%d, NT=%d, Nc=%d, NTH=%d\n', ...
    size(node, 1), size(elem, 1), size(nodeH, 1), size(elemH, 1));

problem = helmholtzLODProblem2D(nodeH, elemH, bdH, node, elem, bdFlag, k, 0, 0);
opts = struct('oversampling', ell, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', useParfor, ...
    'constraintTolerance', 1e-12, 'dropDependentConstraints', true);

t = tic;
P = problem.transfer();
stage(sprintf('transfer P nnz=%d', nnz(P)), t);

t = tic;
Q = problem.interpolation();
stage(sprintf('weighted Clement Q nnz=%d', nnz(Q)), t);

t = tic;
[A, b] = problem.form.global();
stage(sprintf('global form nnz(A)=%d', nnz(A)), t);

t = tic;
patch = lodBuildPatches(nodeH, elemH, node, elem, problem.bdFlagFine, ell, ...
    struct('storeSubmeshes', false));
stage('lightweight patch build', t);

Nf = size(node, 1);
Nc = size(nodeH, 1);
NT = size(elemH, 1);
corrData = cell(NT, 1);
corrStarData = cell(NT, 1);
stats = repmat(emptyStats(), NT, 1);

t = tic;
if useParfor
    parfor T = 1:NT
        [corrData{T}, corrStarData{T}, stats(T)] = ...
            lodElementCorrectorTriplets(T, elemH, problem, opts, P, Q, patch);
    end
else
    for T = 1:NT
        [corrData{T}, corrStarData{T}, stats(T)] = ...
            lodElementCorrectorTriplets(T, elemH, problem, opts, P, Q, patch);
    end
end
stage(sprintf('corrector loop triplets=%d tripletsStar=%d', ...
    sum(cellfun(@tripletCount, corrData)), ...
    sum(cellfun(@tripletCount, corrStarData))), t);
fprintf('Corrector elapsed stats: min %.3g, mean %.3g, max %.3g seconds\n', ...
    min([stats.elapsed]), mean([stats.elapsed]), max([stats.elapsed]));
fprintf('Patch stats: max patchDof %d, max freeDof %d, max constraints %d\n', ...
    max([stats.patchDof]), max([stats.freeDof]), max([stats.constraints]));

t = tic;
trialBasis = lodAssembleCorrectedBasis(P, corrData, Nf, Nc);
clear corrData
testBasis = lodAssembleCorrectedBasis(P, corrStarData, Nf, Nc);
clear corrStarData
stage(sprintf('basis accumulation nnzTrial=%d nnzTest=%d', ...
    nnz(trialBasis), nnz(testBasis)), t);

t = tic;
AH = testBasis' * A * trialBasis;
bH = testBasis' * b;
stage(sprintf('coarse product nnzAH=%d normbH=%.3e', nnz(AH), norm(bH)), t);

fprintf('========== LOD internal stage timing probe complete ==========\n');


function stage(label, timerValue)
fprintf('STAGE %-70s seconds=%10.3f\n', label, toc(timerValue));
printMemory(label);
end


function printMemory(label)
try
    m = memory;
    fprintf('MEM %-72s MATLAB %.3f GB, system available %.3f GB\n', ...
        label, m.MemUsedMATLAB / 2^30, m.PhysicalMemory.Available / 2^30);
catch
end
end


function pool = ensureExactParpool(workers)
workers = max(1, min(round(workers), feature('numcores')));
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', workers);
elseif pool.NumWorkers ~= workers
    delete(pool);
    pool = parpool('local', workers);
end
if pool.NumWorkers ~= workers
    error('debug_lod_stage_timing_probe:workers', ...
        'Requested %d workers, active p.NumWorkers=%d.', workers, pool.NumWorkers);
end
end


function stats = emptyStats()
stats = struct('patchDof', 0, 'freeDof', 0, 'constraints', 0, ...
    'targetDof', 0, 'primalResidual', NaN, 'adjointResidual', NaN, ...
    'constraintResidual', NaN, 'adjointConstraintResidual', NaN, ...
    'elapsed', NaN);
end


function n = tripletCount(data)
if isempty(data)
    n = 0;
else
    n = numel(data.value);
end
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
