function [useParfor, info] = adaptiveParallelWorkerGate(taskName, useParfor, nTasks, opts, preSolveFcn)
% ADAPTIVEPARALLELWORKERGATE  Choose parpool size from one measured subproblem.

if nargin < 1 || isempty(taskName), taskName = 'parallel task'; end
if nargin < 2 || isempty(useParfor), useParfor = false; end
if nargin < 3 || isempty(nTasks), nTasks = 0; end
if nargin < 4, opts = struct(); end
if nargin < 5, preSolveFcn = []; end

opts = normalizeAdaptiveOptions(opts);
useParforRequested = logical(useParfor);
info = emptyGateInfo(taskName, useParforRequested, nTasks, opts);

if ~opts.enabled
    info.notes = appendNote(info.notes, 'adaptive policy disabled');
    return;
end
if ~useParforRequested
    info.notes = appendNote(info.notes, 'parfor not requested');
    return;
end
if nTasks <= 1
    useParfor = false;
    info.useParforEffective = false;
    info.chosenWorkers = 1;
    info.notes = appendNote(info.notes, 'single task uses serial setup');
    return;
end

info.beforeMemory = memorySnapshot();
if isa(preSolveFcn, 'function_handle')
    timer = tic;
    preSolveFcn();
    info.preSolveElapsed = toc(timer);
else
    info.preSolveElapsed = 0;
    info.notes = appendNote(info.notes, 'no representative pre-solve function supplied');
end
info.afterMemory = memorySnapshot();

deltaGB = info.afterMemory.privateCommittedGB - ...
    info.beforeMemory.privateCommittedGB;
if ~isfinite(deltaGB) || deltaGB < 0
    deltaGB = 0;
end
info.preSolveMemoryDeltaGB = deltaGB;

estimatedGB = opts.estimatedWorkerGB;
if ~isfinite(estimatedGB) || estimatedGB < 0
    estimatedGB = 0;
end
info.perWorkerGB = max([opts.perWorkerFloorGB, estimatedGB, ...
    opts.safetyFactor * deltaGB]);

requestedWorkers = requestedWorkerCount(opts, nTasks);
info.requestedWorkers = requestedWorkers;

[memoryWorkers, physicalBudgetGB, commitBudgetGB] = ...
    memoryLimitedWorkerCount(info.afterMemory, info.perWorkerGB, opts);
info.memoryWorkers = memoryWorkers;
info.physicalBudgetGB = physicalBudgetGB;
info.commitBudgetGB = commitBudgetGB;

chosenWorkers = min(requestedWorkers, memoryWorkers);
chosenWorkers = max(1, chosenWorkers);
info.chosenWorkers = chosenWorkers;

if opts.dryRun
    useParfor = false;
    info.useParforEffective = false;
    info.notes = appendNote(info.notes, 'dry run; no pool started');
    return;
end

if chosenWorkers < 2 && opts.forceSerialWhenOneWorker
    useParfor = false;
    info.useParforEffective = false;
    info.notes = appendNote(info.notes, 'memory gate selected one worker; using serial setup');
    return;
end

if opts.startPool
    pool = ensureParallelPool(opts.profile, chosenWorkers);
    info.actualWorkers = pool.NumWorkers;
else
    pool = gcp('nocreate');
    if isempty(pool)
        useParfor = false;
        info.useParforEffective = false;
        info.notes = appendNote(info.notes, ...
            'startPool=false and no existing pool; using serial setup');
        return;
    end
    info.actualWorkers = pool.NumWorkers;
    if pool.NumWorkers ~= chosenWorkers
        info.notes = appendNote(info.notes, ...
            'existing pool size differs from adaptive worker count');
    end
end

useParfor = true;
info.useParforEffective = true;
end


function opts = normalizeAdaptiveOptions(opts)
if islogical(opts) || (isnumeric(opts) && isscalar(opts))
    opts = struct('enabled', logical(opts));
elseif ~isstruct(opts)
    opts = struct();
end

defaults = struct();
defaults.enabled = false;
defaults.profile = 'local';
defaults.startPool = true;
defaults.dryRun = false;
defaults.maxWorkers = [];
defaults.minWorkers = 1;
defaults.safetyFactor = 2;
defaults.clientReserveGB = 8;
defaults.osReserveGB = 16;
defaults.sharedReserveGB = 0;
defaults.outputReserveGB = 0;
defaults.perWorkerFloorGB = 0.25;
defaults.estimatedWorkerGB = NaN;
defaults.forceSerialWhenOneWorker = true;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
opts.enabled = logical(opts.enabled);
opts.startPool = logical(opts.startPool);
opts.dryRun = logical(opts.dryRun);
opts.forceSerialWhenOneWorker = logical(opts.forceSerialWhenOneWorker);
opts.minWorkers = max(1, floor(opts.minWorkers));
if isempty(opts.maxWorkers)
    opts.maxWorkers = feature('numcores');
else
    opts.maxWorkers = max(1, floor(opts.maxWorkers));
end
end


function info = emptyGateInfo(taskName, useParforRequested, nTasks, opts)
info = struct();
info.enabled = opts.enabled;
info.taskName = taskName;
info.nTasks = nTasks;
info.useParforRequested = useParforRequested;
info.useParforEffective = useParforRequested;
info.dryRun = opts.dryRun;
info.requestedWorkers = NaN;
info.memoryWorkers = NaN;
info.chosenWorkers = NaN;
info.actualWorkers = NaN;
info.perWorkerGB = NaN;
info.physicalBudgetGB = NaN;
info.commitBudgetGB = NaN;
info.preSolveElapsed = NaN;
info.preSolveMemoryDeltaGB = NaN;
info.beforeMemory = struct();
info.afterMemory = struct();
info.reserves = struct('clientGB', opts.clientReserveGB, ...
    'osGB', opts.osReserveGB, 'sharedGB', opts.sharedReserveGB, ...
    'outputGB', opts.outputReserveGB, 'safetyFactor', opts.safetyFactor);
info.notes = '';
end


function n = requestedWorkerCount(opts, nTasks)
n = min([opts.maxWorkers, feature('numcores'), nTasks]);
n = max(opts.minWorkers, floor(n));
end


function [nWorkers, physicalBudgetGB, commitBudgetGB] = ...
    memoryLimitedWorkerCount(snapshot, perWorkerGB, opts)
reserveGB = opts.clientReserveGB + opts.osReserveGB + ...
    opts.sharedReserveGB + opts.outputReserveGB;

physicalBudgetGB = snapshot.availablePhysicalGB - reserveGB;
if ~isfinite(physicalBudgetGB)
    physicalBudgetGB = Inf;
end

commitBudgetGB = snapshot.commitHeadroomGB - reserveGB;
if ~isfinite(commitBudgetGB)
    commitBudgetGB = Inf;
end

budgetGB = min(physicalBudgetGB, commitBudgetGB);
if ~isfinite(budgetGB)
    nWorkers = feature('numcores');
elseif budgetGB <= 0
    nWorkers = 1;
else
    nWorkers = floor(budgetGB / max(perWorkerGB, eps));
    nWorkers = max(1, nWorkers);
end
end


function snapshot = memorySnapshot()
snapshot = struct();
snapshot.privateCommittedGB = NaN;
snapshot.availablePhysicalGB = NaN;
snapshot.totalPhysicalGB = NaN;
snapshot.commitHeadroomGB = NaN;
snapshot.commitLimitGB = NaN;
snapshot.maxArrayGB = NaN;

try
    [userMem, systemMem] = memory;
catch
    return;
end

snapshot.privateCommittedGB = bytesToGB(fieldValue(userMem, 'MemUsedMATLAB'));
snapshot.commitHeadroomGB = bytesToGB(fieldValue(userMem, 'MemAvailableAllArrays'));
snapshot.maxArrayGB = bytesToGB(fieldValue(userMem, 'MaxPossibleArrayBytes'));
snapshot.availablePhysicalGB = bytesToGB(nestedFieldValue(systemMem, ...
    {'PhysicalMemory', 'Available'}));
snapshot.totalPhysicalGB = bytesToGB(nestedFieldValue(systemMem, ...
    {'PhysicalMemory', 'Total'}));
if isfinite(snapshot.privateCommittedGB) && isfinite(snapshot.commitHeadroomGB)
    snapshot.commitLimitGB = snapshot.privateCommittedGB + ...
        snapshot.commitHeadroomGB;
end
end


function val = fieldValue(s, name)
val = NaN;
if isstruct(s) && isfield(s, name) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end


function val = nestedFieldValue(s, names)
val = NaN;
for i = 1:numel(names)
    if ~isstruct(s) || ~isfield(s, names{i})
        return;
    end
    s = s.(names{i});
end
if isnumeric(s) && isscalar(s)
    val = double(s);
end
end


function gb = bytesToGB(bytes)
if isfinite(bytes)
    gb = double(bytes) / 2^30;
else
    gb = NaN;
end
end


function pool = ensureParallelPool(profile, nWorkers)
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= nWorkers
    delete(pool);
    pool = [];
end
if isempty(pool)
    if isempty(profile)
        pool = parpool(nWorkers);
    else
        pool = parpool(profile, nWorkers);
    end
end
if pool.NumWorkers ~= nWorkers
    error('adaptiveParallelWorkerGate:poolSize', ...
        'Expected %d workers, but the active pool has %d.', ...
        nWorkers, pool.NumWorkers);
end
end


function txt = appendNote(txt, note)
if isempty(txt)
    txt = note;
else
    txt = sprintf('%s; %s', txt, note);
end
end
