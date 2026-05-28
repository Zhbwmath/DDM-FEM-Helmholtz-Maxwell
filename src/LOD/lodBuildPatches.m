function patch = lodBuildPatches(nodeH, elemH, nodeh, elemh, bdFlagh, ell)
% LODBUILDPATCHES  Aggregate coarse-element LOD patch metadata.

dim = size(nodeH, 2);
if size(nodeh, 2) ~= dim
    error('lodBuildPatches:dim', 'Coarse and fine meshes must have the same dimension.');
end

switch dim
    case 2
        patch = lodBuildPatches2D(nodeH, elemH, nodeh, elemh, bdFlagh, ell);
    case 3
        patch = lodBuildPatches3D(nodeH, elemH, nodeh, elemh, bdFlagh, ell);
    otherwise
        error('lodBuildPatches:dim', 'Only 2D and 3D P1 meshes are supported.');
end
end
