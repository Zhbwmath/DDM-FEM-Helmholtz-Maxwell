function edgeParts = nedelecSubdomainEdges2D(elem, bdFlag, parts)
% NEDELECSUBDOMAINEDGES2D  Convert triangular subdomains to NE_1 edge DOFs.

[edgeElemInc, edge] = nedelecEdgeElementIncidence2D(elem);
[freeEdges, bdEdges] = nedelecFreeEdges2D(elem, bdFlag);

NE = size(edge, 1);
totalIncident = full(sum(edgeElemInc, 2));
freeMask = false(NE, 1);
freeMask(freeEdges) = true;

nSub = numel(parts);
edgeParts = struct();
for s = 1:nSub
    elemIdx = parts(s).elemIdx(:);
    localIncident = full(sum(edgeElemInc(:, elemIdx), 2));
    subEdges = find(localIncident > 0);
    interiorEdges = find(localIncident == totalIncident & freeMask);

    edgeParts(s).elemIdx = elemIdx;
    edgeParts(s).edgeIdx = subEdges;
    edgeParts(s).interiorEdgeIdx = interiorEdges;
    edgeParts(s).boundaryEdgeIdx = setdiff(subEdges, interiorEdges, 'stable');
    edgeParts(s).freeEdgeIdx = intersect(subEdges, freeEdges, 'stable');
end

edgeParts(1).edge = edge;
edgeParts(1).freeEdges = freeEdges;
edgeParts(1).bdEdges = bdEdges;
end
