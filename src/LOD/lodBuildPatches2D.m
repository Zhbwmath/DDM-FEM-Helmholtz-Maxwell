function patch = lodBuildPatches2D(nodeH, elemH, nodeh, elemh, bdFlagh, ell)
% LODBUILDPATCHES2D  Aggregate 2D coarse-element LOD patch metadata.

if nargin < 6 || isempty(ell), ell = 1; end
if ell < 0 || ell ~= floor(ell)
    error('lodBuildPatches2D:ell', 'Oversampling ell must be a nonnegative integer.');
end

NTc = size(elemH, 1);
NTHNodes = size(nodeH, 1);
centroid = (nodeh(elemh(:,1), :) + nodeh(elemh(:,2), :) + nodeh(elemh(:,3), :)) / 3;
[owner, ~] = locateSimplexP1(nodeH, elemH, centroid, 1e-10);
if any(owner == 0)
    bad = find(owner == 0, 1);
    error('lodBuildPatches2D:notNested', ...
        'Fine element centroid %d was not found in the coarse mesh.', bad);
end

fineElemIdsByCoarse = cell(NTc, 1);
for T = 1:NTc
    fineElemIdsByCoarse{T} = find(owner == T);
end

ii = repmat((1:NTc)', 3, 1);
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
    sub = lodSubmesh2D(nodeh, elemh, bdFlagh, fineIds);

    patch.coarseElemIds{T} = coarseIds(:);
    patch.fineElemIds{T} = fineIds(:);
    patch.targetFineElemIds{T} = fineElemIdsByCoarse{T}(:);
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
