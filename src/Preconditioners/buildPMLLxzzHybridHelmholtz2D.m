function precon = buildPMLLxzzHybridHelmholtz2D(node, elem, bdFlag, k, pml, parts, nodeH, elemH, bdH, opts)
% BUILDPMLLXZZHYBRIDHELMHOLTZ2D  Complete PML-LOD LXZZ hybrid instance.

if nargin < 10 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

fine = buildPMLLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, pml, opts.fineOptions);
lod = buildLODHelmholtzPML2D(nodeH, elemH, bdH, node, elem, bdFlag, ...
    k, pml, opts.source, opts.lodOptions);
coarseSpace = buildPMLLODCoarseSpaceHelmholtz2D(fine, lod);
localSolver = buildPMLLxzzLocalSolversHelmholtz2D(fine, parts, opts.localOptions);

hybridOpts = opts.hybridOptions;
hybridOpts.fineSpace = fine;
hybridOpts.coarseSpace = coarseSpace;
hybridOpts.localSolver = localSolver;
hybridOpts.adjointType = opts.adjointType;
hybridOpts.variant = 'pml';

precon = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, hybridOpts);
precon.pml = struct('fineSpace', fine, 'lod', lod, ...
    'localSolver', localSolver, 'coarseSpace', coarseSpace);
end


function opts = localOptions(opts)
defaults = struct();
defaults.source = [];
defaults.fineOptions = struct();
defaults.lodOptions = struct();
defaults.localOptions = struct();
defaults.hybridOptions = struct();
defaults.adjointType = 'reference';
defaults.adaptiveParallelPolicy = false;
defaults.adaptiveParallelOptions = struct();

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
if ~isfield(opts.lodOptions, 'solveCoarse') || isempty(opts.lodOptions.solveCoarse)
    opts.lodOptions.solveCoarse = false;
end
if ~isfield(opts.localOptions, 'localPMLMode') || isempty(opts.localOptions.localPMLMode)
    opts.localOptions.localPMLMode = 'auto';
end
opts.lodOptions = inheritAdaptiveOptions(opts.lodOptions, opts);
opts.localOptions = inheritAdaptiveOptions(opts.localOptions, opts);
opts.hybridOptions = inheritAdaptiveOptions(opts.hybridOptions, opts);
opts.fineOptions.source = opts.source;
end


function child = inheritAdaptiveOptions(child, parent)
if ~isfield(child, 'adaptiveParallelPolicy') || ...
        isempty(child.adaptiveParallelPolicy)
    child.adaptiveParallelPolicy = parent.adaptiveParallelPolicy;
end
if ~isfield(child, 'adaptiveParallelOptions') || ...
        isempty(child.adaptiveParallelOptions)
    child.adaptiveParallelOptions = parent.adaptiveParallelOptions;
end
end
