function local = buildMaxwellDirichletLocalSolvers(fine, parts, opts)
% BUILDMAXWELLDIRICHLETLOCALSOLVERS  Additive Schwarz edge-star Dirichlet solves.

if nargin < 3 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

edgeParts = subdomainEdges(fine, parts);
solverMode = effectiveSolverMode(fine, edgeParts, opts);
nSub = numel(edgeParts);
solvers = cell(nSub, 1);
redIdx = cell(nSub, 1);
extensions = cell(nSub, 1);
solverType = strings(nSub, 1);

if opts.useParfor
    parfor s = 1:nSub
        [solvers{s}, solverType(s), redIdx{s}, extensions{s}] = ...
            setupOneSubdomain(fine, edgeParts(s), solverMode);
    end
else
    for s = 1:nSub
        [solvers{s}, solverType(s), redIdx{s}, extensions{s}] = ...
            setupOneSubdomain(fine, edgeParts(s), solverMode);
    end
end

    function y = applyInverse(r)
        y = zeros(fine.N, size(r, 2));
        if ~isreal(r), y = complex(y); end
        for j = 1:nSub
            if isempty(redIdx{j}), continue; end
            y = y + extensions{j} * solveLocal(solvers{j}, solverType(j), ...
                extensions{j}' * r);
        end
    end

local = struct();
local.applyInverse = @applyInverse;
local.apply = @applyInverse;
local.extensions = extensions;
local.edgeParts = edgeParts;
local.info = localStats(redIdx, solverMode, opts);
end


function opts = localOptions(opts)
defaults = struct();
defaults.solverMode = 'lu';
defaults.useParfor = false;
defaults.localStoredLuLimitGB = 200;
defaults.localLuFillConstant = 40;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function edgeParts = subdomainEdges(fine, parts)
switch fine.dim
    case 2
        edgeParts = nedelecSubdomainEdges2D(fine.elem, fine.bdFlag, parts);
    case 3
        edgeParts = nedelecSubdomainEdges3D(fine.elem, fine.bdFlag, parts);
    otherwise
        error('buildMaxwellDirichletLocalSolvers:dim', ...
            'Only 2D and 3D Maxwell local solvers are supported.');
end
end


function mode = effectiveSolverMode(fine, edgeParts, opts)
mode = lower(opts.solverMode);
if ~strcmp(mode, 'auto')
    return;
end

sizes = zeros(numel(edgeParts), 1);
for s = 1:numel(edgeParts)
    red = reducedInteriorEdges(fine, edgeParts(s));
    sizes(s) = numel(red);
end
bytes = sum(16 * opts.localLuFillConstant .* sizes .* log2(max(sizes, 2)));
estimatedGB = bytes / 1024^3;
if estimatedGB <= opts.localStoredLuLimitGB
    mode = 'lu';
else
    mode = 'direct';
end
end


function redPos = reducedInteriorEdges(fine, edgePart)
interior = edgePart.interiorEdgeIdx(:);
if isempty(interior)
    redPos = zeros(0, 1);
else
    redPos = fine.global2reduced(interior);
    redPos = redPos(redPos > 0);
end
end


function [solver, solverType, redPos, extension] = setupOneSubdomain(fine, edgePart, solverMode)
redPos = reducedInteriorEdges(fine, edgePart);
if isempty(redPos)
    solver = [];
    solverType = "empty";
    extension = sparse(fine.N, 0);
    return;
end

A_loc = fine.A(redPos, redPos);
solverType = string(solverMode);
solver = factorLocalMatrix(A_loc, solverMode);
extension = sparse(redPos, 1:numel(redPos), 1, fine.N, numel(redPos));
end


function solver = factorLocalMatrix(A, solverMode)
switch lower(solverMode)
    case {'direct', 'backslash', 'matrix'}
        solver = struct('A', A);
    case {'lu', 'storedlu'}
        [L, U, p, q] = lu(A, 'vector');
        solver = struct('L', L, 'U', U, 'p', p(:), 'q', q(:));
    otherwise
        error('buildMaxwellDirichletLocalSolvers:solverMode', ...
            'Unknown local solver mode "%s".', solverMode);
end
end


function x = solveLocal(solver, solverType, b)
if solverType == "empty"
    x = zeros(0, size(b, 2));
    return;
end

switch lower(char(solverType))
    case {'direct', 'backslash', 'matrix'}
        x = solver.A \ b;
    case {'lu', 'storedlu'}
        x = zeros(size(b));
        if ~isreal(b), x = complex(x); end
        x(solver.q, :) = solver.U \ (solver.L \ b(solver.p, :));
    otherwise
        error('buildMaxwellDirichletLocalSolvers:solverType', ...
            'Unknown stored solver type "%s".', solverType);
end
end


function info = localStats(redIdx, solverMode, opts)
sizes = cellfun(@numel, redIdx);
if isempty(sizes), sizes = 0; end
info = struct();
info.boundaryCondition = 'Dirichlet artificial boundary on Nedelec edges';
info.nSubdomains = numel(redIdx);
info.localDofMin = min(sizes);
info.localDofMax = max(sizes);
info.localDofMean = mean(sizes);
info.localDofMedian = median(sizes);
info.solverModeRequested = opts.solverMode;
info.solverModeEffective = solverMode;
info.localStoredLuLimitGB = opts.localStoredLuLimitGB;
info.localLuFillConstant = opts.localLuFillConstant;
info.useParfor = opts.useParfor;
end
