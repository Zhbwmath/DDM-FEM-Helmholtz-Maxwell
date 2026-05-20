function [u, convHist] = optimizedSchwarzPoisson2D(node, elem, bdFlag, f, uD, ...
    partitions, alpha, tol, maxIter)
% OPTIMIZEDSCHWARZPOISSON2D  OSM for 2D Poisson (non-overlapping, Robin).
%
%   Robin transmission on interface edges: α ∫_{Γ} u v ds added to stiffness,
%   ∫_{Γ} g v ds added to RHS.
%   partitions must be non-overlapping (overlap=0).

if nargin < 8, alpha = []; end
if nargin < 9, tol = 1e-6; end
if nargin < 10, maxIter = 100; end

nSub = length(partitions);
N = size(node, 1);

if isempty(alpha)
    H = (max(node(:,1)) - min(node(:,1))) / nSub;
    alpha = 0.5 * pi / H;
end

bdNodes = getBoundaryNodes2D(elem, bdFlag);

% Pre-build
A_loc = cell(nSub, 1);  M_loc = cell(nSub, 1);
b_loc = cell(nSub, 1);  u_loc = cell(nSub, 1);
freeLoc = cell(nSub, 1);  bdLoc = cell(nSub, 1);
Mb_ifc = cell(nSub, 1);  % Mb_ifc{s}{ifc} = interface mass matrix

edgeVertPairs = [2 3; 3 1; 1 2];

for s = 1:nSub
    locNode = partitions(s).localNode;
    locElem = partitions(s).localElem;
    nLoc = size(locNode, 1);

    A_loc{s} = assembleStiffness2D(locNode, locElem);
    M_loc{s} = assembleMass2D(locNode, locElem);

    if isnumeric(f)
        b_loc{s} = M_loc{s} * (f * ones(nLoc, 1));
    else
        b_loc{s} = M_loc{s} * f(locNode(:,1), locNode(:,2));
    end

    isBd = false(nLoc, 1);
    isBd(ismember(partitions(s).nodeIdx, bdNodes)) = true;
    bdLoc{s} = find(isBd);
    freeLoc{s} = find(~isBd);

    u_loc{s} = zeros(nLoc, 1);
    if ~isempty(bdLoc{s})
        if isnumeric(uD)
            u_loc{s}(bdLoc{s}) = uD;
        else
            u_loc{s}(bdLoc{s}) = uD(locNode(bdLoc{s},1), locNode(bdLoc{s},2));
        end
    end

    % Assemble interface mass for each interface
    Mb_ifc{s} = cell(partitions(s).nIfaces, 1);
    for ifc = 1:partitions(s).nIfaces
        Mb_ifc{s}{ifc} = assembleEdgeMass(locNode, locElem, ...
            partitions(s).ifaceEdges{ifc}, edgeVertPairs);
    end
end

% Iteration
convHist = zeros(maxIter, 1);

for iter = 1:maxIter
    u_old = u_loc;

    % Compute Robin data from neighbors
    g_iface = cell(nSub, 1);
    for s = 1:nSub
        g_iface{s} = cell(partitions(s).nIfaces, 1);
        for ifc = 1:partitions(s).nIfaces
            sNbr = partitions(s).ifaceNeighbor(ifc);
            nbrIfcNodes = partitions(s).ifaceNodeMap{ifc};
            if isempty(nbrIfcNodes), continue; end

            flux = b_loc{sNbr}(nbrIfcNodes) - ...
                   A_loc{sNbr}(nbrIfcNodes, :) * u_old{sNbr};
            g_iface{s}{ifc} = flux + alpha * u_old{sNbr}(nbrIfcNodes);
        end
    end

    % Solve each subdomain
    for s = 1:nSub
        A_mod = A_loc{s};
        b_mod = b_loc{s};

        for ifc = 1:partitions(s).nIfaces
            if isempty(g_iface{s}{ifc}), continue; end
            Mb = Mb_ifc{s}{ifc};
            gFull = zeros(size(A_mod, 1), 1);
            gFull(partitions(s).ifaceNodes{ifc}) = g_iface{s}{ifc};
            A_mod = A_mod + alpha * Mb;
            b_mod = b_mod + Mb * gFull;
        end

        fL = freeLoc{s};  bL = bdLoc{s};
        if ~isempty(bL)
            u_solve = zeros(size(A_mod, 1), 1);
            if isnumeric(uD)
                u_solve(bL) = uD;
            else
                locNd = partitions(s).localNode;
                u_solve(bL) = uD(locNd(bL,1), locNd(bL,2));
            end
            b_mod(fL) = b_mod(fL) - A_mod(fL, bL) * u_solve(bL);
            u_solve(fL) = A_mod(fL, fL) \ b_mod(fL);
            u_loc{s} = u_solve;
        else
            u_loc{s} = A_mod \ b_mod;
        end
    end

    % Convergence
    jumpNorm = 0;  changeNorm = 0;
    for s = 1:nSub
        for ifc = 1:partitions(s).nIfaces
            sNbr = partitions(s).ifaceNeighbor(ifc);
            if sNbr < s, continue; end
            locNd = partitions(s).ifaceNodes{ifc};
            nbrNd = partitions(s).ifaceNodeMap{ifc};
            if isempty(locNd) || isempty(nbrNd), continue; end
            jumpNorm = max(jumpNorm, norm(u_loc{s}(locNd) - u_loc{sNbr}(nbrNd), inf));
        end
        changeNorm = max(changeNorm, norm(u_loc{s} - u_old{s}, inf));
    end
    convHist(iter) = max(jumpNorm, changeNorm);
    if convHist(iter) < tol
        convHist = convHist(1:iter);  break;
    end
end

u = zeros(N, 1);  count = zeros(N, 1);
for s = 1:nSub
    gIdx = partitions(s).nodeIdx;
    u(gIdx) = u(gIdx) + u_loc{s};
    count(gIdx) = count(gIdx) + 1;
end
mask = count > 0;  u(mask) = u(mask) ./ count(mask);
end


function Mb = assembleEdgeMass(node, elem, ifaceEdges, edgeVertPairs)
% Assemble 1D boundary mass on interface edges (P1).
N = size(node, 1);
nEdges = size(ifaceEdges, 1);
ii = zeros(nEdges * 4, 1);  jj = zeros(nEdges * 4, 1);
ss = zeros(nEdges * 4, 1);  idx = 0;

for i = 1:nEdges
    e = ifaceEdges(i, 1);  k = ifaceEdges(i, 2);
    va = elem(e, edgeVertPairs(k, 1));
    vb = elem(e, edgeVertPairs(k, 2));
    L = norm(node(va, :) - node(vb, :));
    idx = idx + 1; ii(idx) = va; jj(idx) = va; ss(idx) = L/3;
    idx = idx + 1; ii(idx) = vb; jj(idx) = vb; ss(idx) = L/3;
    idx = idx + 1; ii(idx) = va; jj(idx) = vb; ss(idx) = L/6;
    idx = idx + 1; ii(idx) = vb; jj(idx) = va; ss(idx) = L/6;
end
Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
