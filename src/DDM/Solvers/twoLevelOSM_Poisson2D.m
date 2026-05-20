function [u, convHist] = twoLevelOSM_Poisson2D(node, elem, bdFlag, f, uD, ...
    parts, P_H, alpha, tol, maxIter)
% TWOLEVELOSM_POISSON2D  Two-level OSM with coarse correction.
%
%   Each iteration: (1) OSM subdomain solves, (2) coarse correction of residual.
%
%   P_H: prolongation from coarse to fine (on free nodes only).

if nargin < 9, alpha = []; end
if nargin < 10, tol = 1e-6; end
if nargin < 11, maxIter = 100; end

nSub = length(parts);
N = size(node, 1);

if isempty(alpha)
    H = (max(node(:,1)) - min(node(:,1))) / sqrt(nSub);
    alpha = 0.5 * pi / H;
end

bdNodes = getBoundaryNodes2D(elem, bdFlag);

% Assemble global system for residual computation
A_glob = assembleStiffness2D(node, elem);
M_glob = assembleMass2D(node, elem);
if isnumeric(f)
    b_glob = M_glob * (f * ones(N, 1));
else
    b_glob = M_glob * f(node(:,1), node(:,2));
end
freeNodes = setdiff(1:N, bdNodes);
A_ff = A_glob(freeNodes, freeNodes);
b_f = b_glob(freeNodes);
Nf = length(freeNodes);
global2free = zeros(N, 1);
global2free(freeNodes) = (1:Nf)';

% Coarse matrix and solver
R_H = P_H';
A_H = R_H * A_ff * P_H;
try
    L_H = chol(A_H);  coarseSolver = 'chol';
catch
    [L_H, U_H, PLU] = lu(A_H);  coarseSolver = 'lu';
end

% Non-overlapping OSM subdomain setup (reuse existing structure)
edgeVertPairs = [2 3; 3 1; 1 2];
A_loc = cell(nSub, 1);  b_loc = cell(nSub, 1);
u_loc = cell(nSub, 1);  freeLoc = cell(nSub, 1);  bdLoc = cell(nSub, 1);
Mb_ifc = cell(nSub, 1);  locNodeCell = cell(nSub, 1);

for s = 1:nSub
    locNode = parts(s).localNode;
    locElem = parts(s).localElem;
    nLoc = size(locNode, 1);
    locNodeCell{s} = locNode;

    A_loc{s} = assembleStiffness2D(locNode, locElem);
    M_loc_s = assembleMass2D(locNode, locElem);

    if isnumeric(f)
        b_loc{s} = M_loc_s * (f * ones(nLoc, 1));
    else
        b_loc{s} = M_loc_s * f(locNode(:,1), locNode(:,2));
    end

    isBd = false(nLoc, 1);
    isBd(ismember(parts(s).nodeIdx, bdNodes)) = true;
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

    Mb_ifc{s} = cell(parts(s).nIfaces, 1);
    for ifc = 1:parts(s).nIfaces
        Mb_ifc{s}{ifc} = assembleEdgeMass2(locNode, locElem, ...
            parts(s).ifaceEdges{ifc}, edgeVertPairs);
    end
end

% Iteration
convHist = zeros(maxIter, 1);

for iter = 1:maxIter
    u_old = u_loc;

    % ---- Step 1: OSM subdomain solves (parallel Schwarz) -------------------
    g_iface = cell(nSub, 1);
    for s = 1:nSub
        g_iface{s} = cell(parts(s).nIfaces, 1);
        for ifc = 1:parts(s).nIfaces
            sNbr = parts(s).ifaceNeighbor(ifc);
            nbrIfcNodes = parts(s).ifaceNodeMap{ifc};
            if isempty(nbrIfcNodes), continue; end
            flux = b_loc{sNbr}(nbrIfcNodes) - A_loc{sNbr}(nbrIfcNodes, :) * u_old{sNbr};
            g_iface{s}{ifc} = flux + alpha * u_old{sNbr}(nbrIfcNodes);
        end
    end

    for s = 1:nSub
        A_mod = A_loc{s};  b_mod = b_loc{s};
        for ifc = 1:parts(s).nIfaces
            if isempty(g_iface{s}{ifc}), continue; end
            Mb = Mb_ifc{s}{ifc};
            gFull = zeros(size(A_mod, 1), 1);
            gFull(parts(s).ifaceNodes{ifc}) = g_iface{s}{ifc};
            A_mod = A_mod + alpha * Mb;
            b_mod = b_mod + Mb * gFull;
        end
        fL = freeLoc{s};  bL = bdLoc{s};
        if ~isempty(bL)
            u_solve = zeros(size(A_mod, 1), 1);
            if isnumeric(uD)
                u_solve(bL) = uD;
            else
                u_solve(bL) = uD(locNodeCell{s}(bL,1), locNodeCell{s}(bL,2));
            end
            b_mod(fL) = b_mod(fL) - A_mod(fL, bL) * u_solve(bL);
            u_solve(fL) = A_mod(fL, fL) \ b_mod(fL);
            u_loc{s} = u_solve;
        else
            u_loc{s} = A_mod \ b_mod;
        end
    end

    % ---- Step 2: Coarse correction -----------------------------------------
    % Assemble global u from subdomain solutions
    u_glob = zeros(N, 1);  count = zeros(N, 1);
    for s = 1:nSub
        gIdx = parts(s).nodeIdx;
        u_glob(gIdx) = u_glob(gIdx) + u_loc{s};
        count(gIdx) = count(gIdx) + 1;
    end
    u_glob(count > 0) = u_glob(count > 0) ./ count(count > 0);
    u_free = u_glob(freeNodes);

    % Compute residual on free nodes
    r_free = b_f - A_ff * u_free;

    % Coarse correction
    r_H = R_H * r_free;
    if strcmp(coarseSolver, 'chol')
        e_H = L_H \ (L_H' \ r_H);
    else
        e_H = U_H \ (L_H \ (PLU * r_H));
    end
    e_free = P_H * e_H;

    % Update global solution
    u_free = u_free + e_free;
    u_glob = zeros(N, 1);
    u_glob(freeNodes) = u_free;

    % Distribute back to subdomains
    for s = 1:nSub
        gIdx = parts(s).nodeIdx;
        u_loc{s} = u_glob(gIdx);
    end

    % ---- Convergence check -------------------------------------------------
    changeNorm = 0;  jumpNorm = 0;
    for s = 1:nSub
        for ifc = 1:parts(s).nIfaces
            sNbr = parts(s).ifaceNeighbor(ifc);
            if sNbr < s, continue; end
            locNd = parts(s).ifaceNodes{ifc};
            nbrNd = parts(s).ifaceNodeMap{ifc};
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
    gIdx = parts(s).nodeIdx;
    u(gIdx) = u(gIdx) + u_loc{s};
    count(gIdx) = count(gIdx) + 1;
end
mask = count > 0;  u(mask) = u(mask) ./ count(mask);
end


function Mb = assembleEdgeMass2(node, elem, ifaceEdges, edgeVertPairs)
N = size(node, 1);  nEdges = size(ifaceEdges, 1);
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
