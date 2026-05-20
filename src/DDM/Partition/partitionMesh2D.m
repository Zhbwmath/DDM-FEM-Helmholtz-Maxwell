function parts = partitionMesh2D(node, elem, bdFlag, nSub, varargin)
% PARTITIONMESH2D  Partition a 2D triangular mesh into subdomains.
%
%   parts = PARTITIONMESH2D(node, elem, bdFlag, nSub)
%   parts = PARTITIONMESH2D(node, elem, bdFlag, nSub, 'overlap', delta)
%
%   nSub: scalar → strip partition (along x)
%         [nx, ny] → checkerboard grid partition
%
%   Strip:    Ω_s = {e : x_C(e) ∈ strip s}
%   Grid:     Ω_{i,j} = {e : x_C(e) ∈ strip i, y_C(e) ∈ strip j}
%   Overlap extends by delta in each direction.

p = inputParser;
p.addParameter('overlap', 0, @(x) isnumeric(x) && x >= 0);
p.parse(varargin{:});
delta = p.Results.overlap;

% Detect mode
if isscalar(nSub)
    nx = nSub;  ny = 1;  mode = 'strip';
else
    nx = nSub(1);  ny = nSub(2);  mode = 'grid';
end
nTotal = nx * ny;

xMin = min(node(:,1));  xMax = max(node(:,1));
yMin = min(node(:,2));  yMax = max(node(:,2));
NT = size(elem, 1);
Hx = (xMax - xMin) / nx;
Hy = (yMax - yMin) / ny;

% Element centroids
xC = (node(elem(:,1), 1) + node(elem(:,2), 1) + node(elem(:,3), 1)) / 3;
yC = (node(elem(:,1), 2) + node(elem(:,2), 2) + node(elem(:,3), 2)) / 3;

bdNodes = getBoundaryNodes2D(elem, bdFlag);

% ---- Subdomain grid indices to linear index --------------------------------
% sub2lin(i,j) = (j-1)*nx + i
sub2lin = @(i,j) (j-1)*nx + i;

% ---- Non-overlapping base --------------------------------------------------
baseElem = cell(nTotal, 1);
for j = 1:ny
    for i = 1:nx
        s = sub2lin(i, j);
        xL = xMin + (i-1) * Hx;  xR = xMin + i * Hx;
        yL = yMin + (j-1) * Hy;  yR = yMin + j * Hy;
        if i == 1,  xL = xMin - 1e-12; end
        if i == nx, xR = xMax + 1e-12; end
        if j == 1,  yL = yMin - 1e-12; end
        if j == ny, yR = yMax + 1e-12; end
        baseElem{s} = find(xC >= xL & xC < xR & yC >= yL & yC < yR);
    end
end

% ---- Overlapping -----------------------------------------------------------
if delta > 0
    extElem = cell(nTotal, 1);
    for j = 1:ny
        for i = 1:nx
            s = sub2lin(i, j);
            xL = max(xMin, xMin + (i-1) * Hx - delta);
            xR = min(xMax, xMin + i * Hx + delta);
            yL = max(yMin, yMin + (j-1) * Hy - delta);
            yR = min(yMax, yMin + j * Hy + delta);
            extElem{s} = find(xC >= xL & xC <= xR & yC >= yL & yC <= yR);
        end
    end
    useElem = extElem;
else
    useElem = baseElem;
end

% ---- Build global edge list (for fast interior/boundary detection) ----------
edgeVertPairs = [2 3; 3 1; 1 2];
edgesAll = [elem(:, [edgeVertPairs(1,1), edgeVertPairs(1,2)])
            elem(:, [edgeVertPairs(2,1), edgeVertPairs(2,2)])
            elem(:, [edgeVertPairs(3,1), edgeVertPairs(3,2)])];
edgesS = sort(edgesAll, 2);  % sorted (minV, maxV)

% ---- Build partition structs -----------------------------------------------

parts = struct();
for s = 1:nTotal
    eIdx = useElem{s};
    parts(s).elemIdx = eIdx;
    parts(s).nodeIdx = unique(elem(eIdx, :));

    nLocal = length(parts(s).nodeIdx);
    g2l = zeros(size(node, 1), 1);
    g2l(parts(s).nodeIdx) = (1:nLocal)';
    parts(s).localNode = node(parts(s).nodeIdx, :);
    parts(s).localElem = g2l(elem(eIdx, :));
    parts(s).global2local = g2l;

    isBd = ismember(parts(s).nodeIdx, bdNodes);
    parts(s).bdIdx   = find(isBd);
    parts(s).freeIdx = find(~isBd);

    % Interior vs boundary (vectorized via edge detection)
    % Build subdomain edges and count occurrences
    subEdges = [elem(eIdx, [edgeVertPairs(1,1), edgeVertPairs(1,2)])
                elem(eIdx, [edgeVertPairs(2,1), edgeVertPairs(2,2)])
                elem(eIdx, [edgeVertPairs(3,1), edgeVertPairs(3,2)])];
    subEdgesS = sort(subEdges, 2);
    [~, subEid] = ismember(subEdgesS, edgesS, 'rows');
    eCount = accumarray(subEid, 1, [size(edgesS,1), 1]);
    bdPerEdge = eCount(subEid) == 1;  % boundary edges appear once
    % Collect boundary nodes from boundary edges
    bdEdgeNodes = unique(subEdges(bdPerEdge, :));
    bndNodeSet = false(size(node,1), 1);
    bndNodeSet(bdEdgeNodes) = true;
    % Boundary nodes = nodes on boundary edges + global Dirichlet nodes
    bndLocal = find(bndNodeSet(parts(s).nodeIdx));
    intLocal = setdiff((1:nLocal)', bndLocal);
    parts(s).interiorNodeIdx = parts(s).nodeIdx(intLocal);
    parts(s).boundaryNodeIdx = parts(s).nodeIdx(bndLocal);
end

% ---- Interface detection ---------------------------------------------------
for s = 1:nTotal
    parts(s).nIfaces = 0;
    parts(s).ifaceNodes = {};
    parts(s).ifaceEdges = {};
    parts(s).ifaceNeighbor = [];
    parts(s).ifaceNodeMap = {};
end

if delta == 0
    if strcmp(mode, 'strip')
        % Original: check consecutive pairs only
        neighborPairs = [(1:nTotal-1)', (2:nTotal)'];
    else
        % Grid: generate all (i,j) ↔ (i±1,j), (i,j±1) neighbor pairs
        neighborPairs = zeros(0, 2);
        for j = 1:ny
            for i = 1:nx
                s = sub2lin(i, j);
                if i < nx  % right neighbor
                    neighborPairs(end+1, :) = [s, sub2lin(i+1, j)]; %#ok<AGROW>
                end
                if j < ny  % top neighbor
                    neighborPairs(end+1, :) = [s, sub2lin(i, j+1)]; %#ok<AGROW>
                end
            end
        end
    end

    for p = 1:size(neighborPairs, 1)
        sL = neighborPairs(p, 1);  sR = neighborPairs(p, 2);
        sharedNodes = intersect(parts(sL).nodeIdx, parts(sR).nodeIdx);
        if isempty(sharedNodes), continue; end

        for side = 1:2
            if side == 1, sCur = sL; sNbr = sR;
            else,         sCur = sR; sNbr = sL; end

            g2l = parts(sCur).global2local;
            locIfcNodes = g2l(sharedNodes);
            locIfcNodes = locIfcNodes(locIfcNodes > 0);
            locIfcSet = false(size(parts(sCur).localNode, 1), 1);
            locIfcSet(locIfcNodes) = true;

            locElem = parts(sCur).localElem;
            ifaceEdges = zeros(0, 2);
            for ei = 1:size(locElem, 1)
                for k = 1:3
                    va = locElem(ei, edgeVertPairs(k, 1));
                    vb = locElem(ei, edgeVertPairs(k, 2));
                    if locIfcSet(va) && locIfcSet(vb)
                        ifaceEdges(end+1, :) = [ei, k]; %#ok<AGROW>
                    end
                end
            end

            g2lNbr = parts(sNbr).global2local;
            locToNbr = zeros(length(locIfcNodes), 1);
            for i = 1:length(locIfcNodes)
                glb = parts(sCur).nodeIdx(locIfcNodes(i));
                locToNbr(i) = g2lNbr(glb);
            end

            nIf = parts(sCur).nIfaces + 1;
            parts(sCur).nIfaces = nIf;
            parts(sCur).ifaceNodes{nIf} = locIfcNodes;
            parts(sCur).ifaceEdges{nIf} = ifaceEdges;
            parts(sCur).ifaceNeighbor(nIf) = sNbr;
            parts(sCur).ifaceNodeMap{nIf} = locToNbr;
        end
    end
end
end
