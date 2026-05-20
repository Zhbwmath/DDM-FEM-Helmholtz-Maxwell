function [u, convHist] = optimizedSchwarzPoisson1D(node, elem, bdFlag, f, uD, ...
    partitions, alpha, tol, maxIter)
% OPTIMIZEDSCHWARZPOISSON1D  OSM for 1D Poisson (non-overlapping, Robin).
%
%   -u'' = f,  u(0)=u(1)=0
%   Robin transmission: du_i/dn + α u_i = du_j/dn + α u_j  on Γ_{ij}
%
%   partitions must be non-overlapping (overlap=0).

if nargin < 8, alpha = []; end
if nargin < 9, tol = 1e-6; end
if nargin < 10, maxIter = 100; end

nSub = length(partitions);
N = size(node, 1);

if isempty(alpha)
    H = (node(end) - node(1)) / nSub;
    alpha = 0.5 * pi / H;
end

bdNodes = [elem(bdFlag(:,1)==1, 1); elem(bdFlag(:,2)==1, 2)];
bdNodes = unique(bdNodes);

% Pre-build subdomain systems
A_loc = cell(nSub, 1);  M_loc = cell(nSub, 1);
b_loc = cell(nSub, 1);  u_loc = cell(nSub, 1);
freeLoc = cell(nSub, 1);  bdLoc = cell(nSub, 1);

for s = 1:nSub
    locNode = partitions(s).localNode;
    locElem = partitions(s).localElem;
    nLoc = size(locNode, 1);

    A_loc{s} = assembleStiffness1D(locNode, locElem);
    M_loc{s} = assembleMass1D(locNode, locElem);

    if isnumeric(f)
        b_loc{s} = M_loc{s} * (f * ones(nLoc, 1));
    else
        b_loc{s} = M_loc{s} * f(locNode);
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
            u_loc{s}(bdLoc{s}) = uD(locNode(bdLoc{s}));
        end
    end
end

% Iteration
convHist = zeros(maxIter, 1);

for iter = 1:maxIter
    u_old = u_loc;

    % Compute Robin data from neighbors
    gLeft = cell(nSub, 1);  gRight = cell(nSub, 1);
    for s = 1:nSub
        gLeft{s} = 0;  gRight{s} = 0;

        % Left interface: data from left neighbor's right side
        for ifc = 1:partitions(s).nIfaces
            sNbr = partitions(s).ifaceNeighbor(ifc);
            if sNbr > s, continue; end  % left neighbor has index < s
            nbrLocIfc = partitions(s).ifaceNodeMap{ifc};
            if isempty(nbrLocIfc), continue; end
            k = nbrLocIfc(1);  % 1D: single interface node
            flux = b_loc{sNbr}(k) - A_loc{sNbr}(k, :) * u_old{sNbr};
            gLeft{s} = flux + alpha * u_old{sNbr}(k);
        end

        % Right interface: data from right neighbor's left side
        for ifc = 1:partitions(s).nIfaces
            sNbr = partitions(s).ifaceNeighbor(ifc);
            if sNbr < s, continue; end  % right neighbor has index > s
            nbrLocIfc = partitions(s).ifaceNodeMap{ifc};
            if isempty(nbrLocIfc), continue; end
            k = nbrLocIfc(1);
            flux = b_loc{sNbr}(k) - A_loc{sNbr}(k, :) * u_old{sNbr};
            gRight{s} = flux + alpha * u_old{sNbr}(k);
        end
    end

    % Solve each subdomain with Robin BC
    for s = 1:nSub
        A_mod = A_loc{s};
        b_mod = b_loc{s};

        % Add Robin terms at interfaces
        for ifc = 1:partitions(s).nIfaces
            ifcNodes = partitions(s).ifaceNodes{ifc};
            if isempty(ifcNodes), continue; end
            k = ifcNodes(1);  % 1D: single node
            sNbr = partitions(s).ifaceNeighbor(ifc);

            A_mod(k, k) = A_mod(k, k) + alpha;
            if sNbr < s
                b_mod(k) = b_mod(k) + gLeft{s};
            else
                b_mod(k) = b_mod(k) + gRight{s};
            end
        end

        % Dirichlet BC and solve
        fL = freeLoc{s};  bL = bdLoc{s};
        if ~isempty(bL)
            u_solve = zeros(size(A_mod, 1), 1);
            if isnumeric(uD)
                u_solve(bL) = uD;
            else
                u_solve(bL) = uD(partitions(s).localNode(bL));
            end
            b_mod(fL) = b_mod(fL) - A_mod(fL, bL) * u_solve(bL);
            u_solve(fL) = A_mod(fL, fL) \ b_mod(fL);
            u_loc{s} = u_solve;
        else
            u_loc{s} = A_mod \ b_mod;
        end
    end

    % Convergence check
    jumpNorm = 0;  changeNorm = 0;
    for s = 1:nSub
        for ifc = 1:partitions(s).nIfaces
            sNbr = partitions(s).ifaceNeighbor(ifc);
            if sNbr < s, continue; end
            locNd = partitions(s).ifaceNodes{ifc};
            nbrNd = partitions(s).ifaceNodeMap{ifc};
            if isempty(locNd) || isempty(nbrNd), continue; end
            jumpNorm = max(jumpNorm, abs(u_loc{s}(locNd(1)) - u_loc{sNbr}(nbrNd(1))));
        end
        changeNorm = max(changeNorm, norm(u_loc{s} - u_old{s}, inf));
    end
    convHist(iter) = max(jumpNorm, changeNorm);
    if convHist(iter) < tol
        convHist = convHist(1:iter);  break;
    end
end

% Assemble global solution (average at shared nodes)
u = zeros(N, 1);  count = zeros(N, 1);
for s = 1:nSub
    gIdx = partitions(s).nodeIdx;
    u(gIdx) = u(gIdx) + u_loc{s};
    count(gIdx) = count(gIdx) + 1;
end
mask = count > 0;  u(mask) = u(mask) ./ count(mask);
end
