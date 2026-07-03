function sub = lodGetPatchSubmesh(patch, T)
% LODGETPATCHSUBMESH  Return a stored or lazily materialized patch submesh.

if isfield(patch, 'local2global') && numel(patch.local2global) >= T && ...
        ~isempty(patch.local2global{T})
    sub = struct();
    sub.local2global = patch.local2global{T};
    sub.localElem = patch.localElem{T};
    sub.localNode = patch.localNode{T};
    sub.localBdFlag = patch.localBdFlag{T};
    sub.freeLocalDof = patch.freeLocalDof{T};
    sub.boundaryLocalDof = patch.boundaryLocalDof{T};
    sub.artificialBoundaryLocalDof = patch.artificialBoundaryLocalDof{T};
    sub.physicalBoundaryLocalDof = patch.physicalBoundaryLocalDof{T};
    return;
end

if ~isfield(patch, 'node') || ~isfield(patch, 'elem') || ...
        ~isfield(patch, 'bdFlag') || ~isfield(patch, 'fineElemIds')
    error('lodGetPatchSubmesh:missingData', ...
        'Patch does not store submeshes or enough data to build them lazily.');
end

dim = size(patch.node, 2);
switch dim
    case 2
        sub = lodSubmesh2D(patch.node, patch.elem, patch.bdFlag, ...
            patch.fineElemIds{T});
    case 3
        sub = lodSubmesh3D(patch.node, patch.elem, patch.bdFlag, ...
            patch.fineElemIds{T});
    otherwise
        error('lodGetPatchSubmesh:dim', ...
            'Only 2D and 3D P1 patch submeshes are supported.');
end
end
