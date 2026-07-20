function [coarse, info, lod] = buildLODMaxwellCoarseSpace(fine, nodeH, elemH, bdH, opts)
% BUILDLODMAXWELLCOARSESPACE  Normalize reference LOD4Maxwell coarse basis.

if nargin < 5 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

if ~isempty(opts.lod)
    lod = opts.lod;
    referenceInfo = struct('mode', 'injected lod object', ...
        'path', '', 'addedPath', false);
else
    [srcPath, referenceRoot] = resolveReferencePath(opts.lodReferencePath);
    oldPath = path;
    cleanupPath = onCleanup(@() path(oldPath));
    addpath(genpath(srcPath), '-end');
    builderName = sprintf('buildLODMaxwell%dD', fine.dim);
    if exist(builderName, 'file') ~= 2
        error('buildLODMaxwellCoarseSpace:referenceBuilder', ...
            'Could not find %s after adding reference path "%s".', ...
            builderName, srcPath);
    end
    lodOpts = lodOptions(opts.lodOptions);
    builder = str2func(builderName);
    lod = builder(nodeH, elemH, bdH, fine.node, fine.elem, fine.bdFlag, ...
        fine.kappa, [], lodOpts);
    referenceInfo = struct('mode', 'LOD4Maxwell reference path', ...
        'path', referenceRoot, 'srcPath', srcPath, 'addedPath', true);
end

validateLodData(fine, lod, opts);
[part, selectedBasis] = selectLodPart(lod, opts.lodBasis);
[trial, basisSource] = correctedBasisFromLod(lod, part, selectedBasis);
if size(trial, 1) ~= fine.N
    error('buildLODMaxwellCoarseSpace:basisSize', ...
        'LOD basis has %d rows, but the fine Maxwell space has %d free edges.', ...
        size(trial, 1), fine.N);
end

if opts.recomputeCoarseMatrix || ~isfield(part, 'unified') || ...
        ~isfield(part.unified, 'A') || isempty(part.unified.A)
    AH = trial' * fine.A * trial;
    matrixSource = 'recomputed from fine.A';
else
    AH = part.unified.A;
    matrixSource = sprintf('LOD %s unified.A', selectedBasis);
end

coarse = struct();
coarse.trial = trial;
coarse.test = trial;
coarse.nativeTrial = trial;
coarse.nativeTest = trial;
coarse.embedding = speye(fine.N);
coarse.AH = AH;
coarse.solve = @(r) AH \ r;
coarse.solveAdjoint = @(r) euclideanTranspose(AH) \ r;
coarse.energyAdjointTrial = [];
coarse.info = struct('description', ...
    sprintf('Maxwell %s LOD coarse space from reference LOD4Maxwell', selectedBasis), ...
    'sourcePolicy', 'reference-only; no LOD4Maxwell files are committed here', ...
    'reference', referenceInfo, 'selectedBasis', selectedBasis, ...
    'basisSource', basisSource, 'matrixSource', matrixSource);

info = coarse.info;
info.lodDiagnostics = getFieldOrEmpty(lod, 'diagnostics');
info.nCoarseDofs = size(trial, 2);
info.nFineDofs = fine.N;
end


function opts = localOptions(opts)
defaults = struct();
defaults.lod = [];
defaults.lodReferencePath = '';
defaults.lodOptions = struct();
defaults.lodBasis = 'localized';
defaults.recomputeCoarseMatrix = false;
defaults.checkDofConsistency = true;
defaults.checkFineMatrixConsistency = false;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function lodOpts = lodOptions(lodOpts)
if ~isfield(lodOpts, 'solveFine') || isempty(lodOpts.solveFine)
    lodOpts.solveFine = false;
end
if ~isfield(lodOpts, 'buildLocalized') || isempty(lodOpts.buildLocalized)
    lodOpts.buildLocalized = true;
end
if ~isfield(lodOpts, 'storeLocalizedCorrectors') || isempty(lodOpts.storeLocalizedCorrectors)
    lodOpts.storeLocalizedCorrectors = true;
end
if ~isfield(lodOpts, 'storeCorrectedBasis') || isempty(lodOpts.storeCorrectedBasis)
    lodOpts.storeCorrectedBasis = false;
end
if ~isfield(lodOpts, 'localSolverMode') || isempty(lodOpts.localSolverMode)
    lodOpts.localSolverMode = 'direct';
end
end


function [srcPath, referenceRoot] = resolveReferencePath(inputPath)
if isempty(inputPath)
    thisFile = mfilename('fullpath');
    repoRoot = fileparts(fileparts(fileparts(thisFile)));
    referenceRoot = fullfile(fileparts(repoRoot), 'LOD4Maxwell');
else
    referenceRoot = char(inputPath);
end

if strcmpi(getLastPathPart(referenceRoot), 'src')
    srcPath = referenceRoot;
    referenceRoot = fileparts(referenceRoot);
else
    srcPath = fullfile(referenceRoot, 'src');
end

if exist(srcPath, 'dir') ~= 7
    error('buildLODMaxwellCoarseSpace:referencePath', ...
        ['LOD4Maxwell reference src path "%s" does not exist. ', ...
        'Pass opts.lodReferencePath or opts.lod with an already-built object.'], ...
        srcPath);
end
end


function name = getLastPathPart(p)
[~, name] = fileparts(p);
end


function validateLodData(fine, lod, opts)
required = {'data'};
for i = 1:numel(required)
    if ~isfield(lod, required{i})
        error('buildLODMaxwellCoarseSpace:lodObject', ...
            'LOD object is missing field "%s".', required{i});
    end
end
if ~isfield(lod.data, 'A') || ~isfield(lod.data, 'freeEdgeh')
    error('buildLODMaxwellCoarseSpace:lodData', ...
        'LOD data must define A and freeEdgeh.');
end
if size(lod.data.A, 1) ~= fine.N
    error('buildLODMaxwellCoarseSpace:fineSize', ...
        'LOD fine matrix has %d rows, but fine.N is %d.', ...
        size(lod.data.A, 1), fine.N);
end
if opts.checkDofConsistency && ~isequal(lod.data.freeEdgeh(:), fine.freeEdges(:))
    error('buildLODMaxwellCoarseSpace:dofConvention', ...
        'LOD4Maxwell and current repo free-edge sets do not match.');
end
if opts.checkFineMatrixConsistency
    relA = norm(lod.data.A - fine.A, 'fro') / max(1, norm(fine.A, 'fro'));
    if relA > 1e-10
        error('buildLODMaxwellCoarseSpace:fineMatrixMismatch', ...
            'LOD fine matrix differs from current fine.A by %.3e.', relA);
    end
end
end


function [part, selectedBasis] = selectLodPart(lod, requestedBasis)
requested = lower(requestedBasis);
part = [];
selectedBasis = '';

if any(strcmp(requested, {'localized', 'local', 'auto'}))
    if isfield(lod, 'localized') && ~isempty(lod.localized)
        part = lod.localized;
        selectedBasis = 'localized';
    end
end
if isempty(part) && any(strcmp(requested, {'global', 'auto', 'localized', 'local'}))
    if isfield(lod, 'global') && ~isempty(lod.global) && ...
            isfield(lod.global, 'built') && lod.global.built
        part = lod.global;
        selectedBasis = 'global';
    end
end
if isempty(part)
    error('buildLODMaxwellCoarseSpace:lodBasis', ...
        'LOD object does not contain a usable "%s" corrected coarse basis.', ...
        requestedBasis);
end
end


function [basis, source] = correctedBasisFromLod(lod, part, selectedBasis)
if isfield(part, 'unified') && isfield(part.unified, 'basis') && ...
        ~isempty(part.unified.basis)
    basis = part.unified.basis;
    source = sprintf('%s.unified.basis', selectedBasis);
    return;
end

switch selectedBasis
    case 'localized'
        requireField(part, 'X', 'localized basis');
        requireField(part.X, 'corrector', 'localized X corrector');
        requireField(part, 'scalar', 'localized basis');
        requireField(part.scalar, 'correctorPotential', ...
            'localized scalar corrector potential');
        basis = lod.data.Pcurl - part.X.corrector - ...
            lod.data.G * part.scalar.correctorPotential;
        source = 'reconstructed from localized Pcurl, X.corrector, and scalar.correctorPotential';
    case 'global'
        requireField(part, 'X', 'global basis');
        requireField(part.X, 'corrector', 'global X corrector');
        requireField(part, 'scalar', 'global basis');
        if isfield(part.scalar, 'corrector') && ~isempty(part.scalar.corrector)
            gradCorrector = part.scalar.corrector;
        else
            requireField(part.scalar, 'correctorPotential', ...
                'global scalar corrector potential');
            gradCorrector = lod.data.G * part.scalar.correctorPotential;
        end
        basis = lod.data.Pcurl - part.X.corrector - gradCorrector;
        source = 'reconstructed from global Pcurl, X.corrector, and scalar corrector';
    otherwise
        error('buildLODMaxwellCoarseSpace:lodBasis', ...
            'Unsupported selected basis "%s".', selectedBasis);
end
end


function requireField(s, fieldName, context)
if ~isfield(s, fieldName) || isempty(s.(fieldName))
    error('buildLODMaxwellCoarseSpace:missingField', ...
        'Missing %s field "%s".', context, fieldName);
end
end


function value = getFieldOrEmpty(s, fieldName)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = [];
end
end
