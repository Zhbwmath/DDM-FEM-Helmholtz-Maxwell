function patch = lodBuildPatches3D(nodeH, elemH, nodeh, elemh, bdFlagh, ell, opts)
% LODBUILDPATCHES3D  Aggregate 3D coarse-element LOD patch metadata.

if nargin < 6 || isempty(ell), ell = 1; end
if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'storeSubmeshes') || isempty(opts.storeSubmeshes)
    opts.storeSubmeshes = true;
end
if ell < 0 || ell ~= floor(ell)
    error('lodBuildPatches3D:ell', 'Oversampling ell must be a nonnegative integer.');
end

NTc = size(elemH, 1);
NTHNodes = size(nodeH, 1);
centroid = (nodeh(elemh(:,1), :) + nodeh(elemh(:,2), :) + ...
    nodeh(elemh(:,3), :) + nodeh(elemh(:,4), :)) / 4;
[owner, ~] = locateSimplexP1(nodeH, elemH, centroid, 1e-10);
if any(owner == 0)
    bad = find(owner == 0, 1);
    error('lodBuildPatches3D:notNested', ...
        'Fine element centroid %d was not found in the coarse mesh.', bad);
end

fineElemIdsByCoarse = groupIndicesByOwner(owner, NTc);

ii = repmat((1:NTc)', 4, 1);
jj = elemH(:);
T2N = sparse(ii, jj, 1, NTc, NTHNodes);
T2T = spones(T2N * T2N');
if ell == 0
    patchPattern = speye(NTc);
else
    patchPattern = spones((T2T + speye(NTc))^ell);
end

patch = struct();
patch.oversampling = ell;
patch.dimension = 3;
patch.node = nodeh;
patch.elem = elemh;
patch.bdFlag = bdFlagh;
patch.storeSubmeshes = opts.storeSubmeshes;
patch.fineElemIdsByCoarse = fineElemIdsByCoarse;
patch.coarseElemIds = cell(NTc, 1);
patch.fineElemIds = cell(NTc, 1);
patch.targetFineElemIds = cell(NTc, 1);
patch.local2global = cell(NTc, 1);
patch.localElem = cell(NTc, 1);
patch.localNode = cell(NTc, 1);
patch.localBdFlag = cell(NTc, 1);
patch.freeLocalDof = cell(NTc, 1);
patch.boundaryLocalDof = cell(NTc, 1);
patch.artificialBoundaryLocalDof = cell(NTc, 1);
patch.physicalBoundaryLocalDof = cell(NTc, 1);

for T = 1:NTc
    coarseIds = find(patchPattern(T, :));
    fineIds = unique(vertcat(fineElemIdsByCoarse{coarseIds}));

    patch.coarseElemIds{T} = coarseIds(:);
    patch.fineElemIds{T} = fineIds(:);
    patch.targetFineElemIds{T} = fineElemIdsByCoarse{T}(:);
    if opts.storeSubmeshes
        sub = lodSubmesh3D(nodeh, elemh, bdFlagh, fineIds);
        patch.local2global{T} = sub.local2global;
        patch.localElem{T} = sub.localElem;
        patch.localNode{T} = sub.localNode;
        patch.localBdFlag{T} = sub.localBdFlag;
        patch.freeLocalDof{T} = sub.freeLocalDof;
        patch.boundaryLocalDof{T} = sub.boundaryLocalDof;
        patch.artificialBoundaryLocalDof{T} = sub.artificialBoundaryLocalDof;
        patch.physicalBoundaryLocalDof{T} = sub.physicalBoundaryLocalDof;
    end
end
end


function groups = groupIndicesByOwner(owner, nOwner)
owner = owner(:);
[ownerSorted, perm] = sort(owner);
counts = accumarray(ownerSorted, 1, [nOwner, 1]);
offset = [0; cumsum(counts)];
groups = cell(nOwner, 1);
for k = 1:nOwner
    groups{k} = perm(offset(k) + (1:counts(k)));
end
end
