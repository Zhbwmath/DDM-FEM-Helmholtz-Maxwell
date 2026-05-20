function parts = partitionMesh1D(node, elem, bdFlag, nSub, varargin)
% PARTITIONMESH1D  Partition a 1D mesh into subdomains.
%
%   parts = PARTITIONMESH1D(node, elem, bdFlag, nSub)
%   parts = PARTITIONMESH1D(node, elem, bdFlag, nSub, 'overlap', delta)
%
%   Non-overlapping base: split [a,b] into nSub intervals of equal length H.
%     Element e ∈ Ω_s^0 if its centroid lies in [a+(s-1)*H, a+s*H].
%   Overlap (optional): extend each Ω_s^0 by geometric distance delta.
%     Element e ∈ Ω_s if centroid ∈ [a+(s-1)*H-delta, a+s*H+delta] ∩ [a,b].
%
%   Output struct fields:
%     Common:  .elemIdx, .nodeIdx, .localNode, .localElem, .global2local
%     Overlap: .interiorNodeIdx  (nodes strictly inside Ω_i)
%              .boundaryNodeIdx  (nodes on ∂Ω_i, Dirichlet u=0 for ASM)
%     Non-ov:  .freeIdx, .bdIdx  (global Dirichlet boundary)
%              .ifaceNodes, .ifaceEdges, .ifaceNeighbor, .ifaceNodeMap

p = inputParser;
p.addParameter('overlap', 0, @(x) isnumeric(x) && x >= 0);
p.parse(varargin{:});
delta = p.Results.overlap;

a = node(1);  b = node(end);
NT = size(elem, 1);
H = (b - a) / nSub;

% Element centroids
xC = (node(elem(:,1)) + node(elem(:,2))) / 2;

% Global Dirichlet boundary nodes
bdNodes = [elem(bdFlag(:,1)==1, 1); elem(bdFlag(:,2)==1, 2)];
bdNodes = unique(bdNodes);

% ---- Non-overlapping base --------------------------------------------------
baseElem = cell(nSub, 1);
for s = 1:nSub
    xLeft  = a + (s-1) * H;
    xRight = a + s * H;
    if s == 1,     xLeft = a - 1e-12; end
    if s == nSub,  xRight = b + 1e-12; end
    baseElem{s} = find(xC >= xLeft & xC < xRight);
end

% ---- Overlapping (if delta > 0) -------------------------------------------
if delta > 0
    extElem = cell(nSub, 1);
    for s = 1:nSub
        xLeft  = max(a, a + (s-1) * H - delta);
        xRight = min(b, a + s * H + delta);
        extElem{s} = find(xC >= xLeft & xC <= xRight);
    end
    useElem = extElem;
else
    useElem = baseElem;
end

% ---- Build partition structs -----------------------------------------------
parts = struct();
for s = 1:nSub
    eIdx = useElem{s};
    parts(s).elemIdx = eIdx;
    parts(s).nodeIdx = unique(elem(eIdx, :));

    nLocal = length(parts(s).nodeIdx);
    g2l = zeros(size(node, 1), 1);
    g2l(parts(s).nodeIdx) = (1:nLocal)';
    parts(s).localNode = node(parts(s).nodeIdx);
    parts(s).localElem = g2l(elem(eIdx, :));
    parts(s).global2local = g2l;

    % Classify global Dirichlet boundary nodes
    isBd = ismember(parts(s).nodeIdx, bdNodes);
    parts(s).bdIdx   = find(isBd);
    parts(s).freeIdx = find(~isBd);

    % Interior vs boundary nodes for ASM (always computed)
    % Interior: all elements containing this node are in Ω_i
    % Boundary: at least one element containing this node is outside Ω_i
    subNodes = parts(s).nodeIdx;
    interior = [];  boundary = [];
    for ni = 1:length(subNodes)
        gNode = subNodes(ni);
        elemWithNode = find(elem(:,1) == gNode | elem(:,2) == gNode);
        outsideSub = setdiff(elemWithNode, eIdx);
        if isempty(outsideSub)
            interior = [interior; gNode]; %#ok<AGROW>
        else
            boundary = [boundary; gNode]; %#ok<AGROW>
        end
    end
    parts(s).interiorNodeIdx = interior;
    parts(s).boundaryNodeIdx = boundary;
end

% ---- Interface detection (non-overlapping) ---------------------------------
for s = 1:nSub
    parts(s).nIfaces = 0;
    parts(s).ifaceNodes = {};
    parts(s).ifaceEdges = {};
    parts(s).ifaceNeighbor = [];
    parts(s).ifaceNodeMap = [];
end

if delta == 0
    for s = 1:nSub-1
        sharedNodes = intersect(parts(s).nodeIdx, parts(s+1).nodeIdx);
        if isempty(sharedNodes), continue; end

        for side = 1:2
            if side == 1, sCur = s; sNbr = s+1;
            else,         sCur = s+1; sNbr = s; end

            g2l = parts(sCur).global2local;
            locIfc = g2l(sharedNodes);
            locIfc = locIfc(locIfc > 0);

            g2lNbr = parts(sNbr).global2local;
            locToNbr = arrayfun(@(g) g2lNbr(g), parts(sCur).nodeIdx(locIfc));

            nIf = parts(sCur).nIfaces + 1;
            parts(sCur).nIfaces = nIf;
            parts(sCur).ifaceNodes{nIf} = locIfc;
            parts(sCur).ifaceEdges{nIf} = [];
            parts(sCur).ifaceNeighbor(nIf) = sNbr;
            parts(sCur).ifaceNodeMap{nIf} = locToNbr;
        end
    end
end
end
