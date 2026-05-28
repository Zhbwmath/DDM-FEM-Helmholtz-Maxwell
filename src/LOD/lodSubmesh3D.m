function sub = lodSubmesh3D(node, elem, bdFlag, elemIds)
% LODSUBMESH3D  Build local tetrahedral patch mesh and boundary DOF sets.

local2global = unique(elem(elemIds, :));
global2local = zeros(size(node, 1), 1);
global2local(local2global) = 1:numel(local2global);

localElem = global2local(elem(elemIds, :));
localNode = node(local2global, :);
localBdFlag = bdFlag(elemIds, :);

faceVertex = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
[boundaryLocalDof, boundaryFaces] = localBoundaryFaces(localElem, faceVertex);

physicalBoundaryLocalDof = [];
physicalBoundaryFaces = zeros(0, 3);
for k = 1:4
    isBd = localBdFlag(:, k) == 1;
    if any(isBd)
        facesK = localElem(isBd, faceVertex(k, :));
        physicalBoundaryFaces = [physicalBoundaryFaces; facesK]; %#ok<AGROW>
        physicalBoundaryLocalDof = [physicalBoundaryLocalDof; facesK]; %#ok<AGROW>
    end
end
physicalBoundaryLocalDof = unique(physicalBoundaryLocalDof(:));

if isempty(physicalBoundaryFaces)
    isPhysicalBoundaryFace = false(size(boundaryFaces, 1), 1);
else
    isPhysicalBoundaryFace = ismember(sort(boundaryFaces, 2), ...
        sort(physicalBoundaryFaces, 2), 'rows');
end
artificialBoundaryFaces = boundaryFaces(~isPhysicalBoundaryFace, :);
artificialBoundaryLocalDof = unique(artificialBoundaryFaces(:));
isFree = true(numel(local2global), 1);
isFree(artificialBoundaryLocalDof) = false;

sub = struct();
sub.local2global = local2global(:);
sub.localElem = localElem;
sub.localNode = localNode;
sub.localBdFlag = localBdFlag;
sub.freeLocalDof = isFree;
sub.boundaryLocalDof = boundaryLocalDof(:);
sub.artificialBoundaryLocalDof = artificialBoundaryLocalDof(:);
sub.physicalBoundaryLocalDof = physicalBoundaryLocalDof(:);
sub.boundaryFace = boundaryFaces;
end


function [boundaryLocalDof, boundaryFaces] = localBoundaryFaces(localElem, faceVertex)
NT = size(localElem, 1);
allFaces = zeros(4 * NT, 3);
for k = 1:4
    rows = (k - 1) * NT + (1:NT);
    allFaces(rows, :) = localElem(:, faceVertex(k, :));
end

sortedFaces = sort(allFaces, 2);
[~, ~, faceId] = unique(sortedFaces, 'rows');
counts = accumarray(faceId, 1);
isBoundary = counts(faceId) == 1;
boundaryFaces = allFaces(isBoundary, :);
boundaryLocalDof = unique(boundaryFaces(:));
end
