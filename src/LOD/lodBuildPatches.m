function patch = lodBuildPatches(nodeH, elemH, nodeh, elemh, bdFlagh, ell, opts)
% LODBUILDPATCHES  Aggregate coarse-element LOD patch metadata.

if nargin < 7 || isempty(opts), opts = struct(); end
dim = size(nodeH, 2);
if size(nodeh, 2) ~= dim
    error('lodBuildPatches:dim', 'Coarse and fine meshes must have the same dimension.');
end

switch dim
    case 2
        patch = lodBuildPatches2D(nodeH, elemH, nodeh, elemh, bdFlagh, ell, opts);
    case 3
        patch = lodBuildPatches3D(nodeH, elemH, nodeh, elemh, bdFlagh, ell, opts);
    otherwise
        error('lodBuildPatches:dim', 'Only 2D and 3D P1 meshes are supported.');
end
end
