function sub = lodSubmesh2D(node, elem, bdFlag, elemIds)
% LODSUBMESH2D  Build local patch mesh and artificial-boundary DOF sets.

local2global = unique(elem(elemIds, :));
global2local = zeros(size(node, 1), 1);
global2local(local2global) = 1:numel(local2global);

localElem = global2local(elem(elemIds, :));
localNode = node(local2global, :);
localBdFlag = bdFlag(elemIds, :);

[edge, ~, ~, edge2elem] = edgeMesh2D(localElem);
boundaryEdges = edge(edge2elem(:,2) == 0, :);
boundaryLocalDof = unique(boundaryEdges(:));

edgeVertex = [2 3; 3 1; 1 2];
physicalBoundaryLocalDof = [];
physicalBoundaryEdges = zeros(0, 2);
for k = 1:3
    isBd = localBdFlag(:, k) == 1;
    if any(isBd)
        edgesK = localElem(isBd, edgeVertex(k, :));
        physicalBoundaryEdges = [physicalBoundaryEdges; edgesK]; %#ok<AGROW>
        physicalBoundaryLocalDof = [physicalBoundaryLocalDof; edgesK]; %#ok<AGROW>
    end
end
physicalBoundaryLocalDof = unique(physicalBoundaryLocalDof(:));

if isempty(physicalBoundaryEdges)
    isPhysicalBoundaryEdge = false(size(boundaryEdges, 1), 1);
else
    isPhysicalBoundaryEdge = ismember(sort(boundaryEdges, 2), ...
        sort(physicalBoundaryEdges, 2), 'rows');
end
artificialBoundaryEdges = boundaryEdges(~isPhysicalBoundaryEdge, :);
artificialBoundaryLocalDof = unique(artificialBoundaryEdges(:));
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
end
