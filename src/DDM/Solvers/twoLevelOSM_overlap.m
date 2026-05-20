function [u, convHist] = twoLevelOSM_overlap(node, elem, bdFlag, f, uD, ...
    parts, P_H, alpha, tol, maxIter)
% TWOLEVELOSM_OVERLAP  Two-level OSM with overlapping subdomains.
%
%   Each iteration: (1) overlapping OSM, then (2) coarse correction.
%   Same as twoLevelOSM_Poisson2D but uses overlapping OSM for fine level.

if nargin < 9, alpha = []; end
if nargin < 10, tol = 1e-6; end
if nargin < 11, maxIter = 100; end

nSub = length(parts);
N = size(node, 1);
NT = size(elem, 1);

if isempty(alpha)
    H = (max(node(:,1)) - min(node(:,1))) / sqrt(nSub);
    alpha = pi / H;
end

bdNodes = getBoundaryNodes2D(elem, bdFlag);

% Global system for residual computation
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

% Coarse solver
R_H = P_H';
A_H = R_H * A_ff * P_H;
try
    L_H = chol(A_H);  coarseSolver = 'chol';
catch
    [L_H, U_H, PLU] = lu(A_H);  coarseSolver = 'lu';
end

% ---- Precompute interface edges (same as overlapping OSM) -----------------
edgeVP = [2 3; 3 1; 1 2];
edge2elem = containers.Map('KeyType', 'char', 'ValueType', 'any');
for e = 1:NT
    for k = 1:3
        va = elem(e, edgeVP(k, 1));
        vb = elem(e, edgeVP(k, 2));
        key = sprintf('%d,%d', min(va,vb), max(va,vb));
        if isKey(edge2elem, key)
            edge2elem(key) = [edge2elem(key), e];
        else
            edge2elem(key) = e;
        end
    end
end

edgeInfo = cell(nSub, 1);
for s = 1:nSub
    g2l = parts(s).global2local;
    eIdx = parts(s).elemIdx;
    eSet = false(NT, 1);  eSet(eIdx) = true;
    locElem = parts(s).localElem;

    eInfo = struct('locEdgeNodes', {}, 'elemOut', {}, 'normal', {}, 'length', {});
    nEdges = 0;

    for ei = 1:size(locElem, 1)
        eGlob = eIdx(ei);
        for k = 1:3
            vaGlob = elem(eGlob, edgeVP(k, 1));
            vbGlob = elem(eGlob, edgeVP(k, 2));
            key = sprintf('%d,%d', min(vaGlob, vbGlob), max(vaGlob, vbGlob));
            adjElems = edge2elem(key);
            otherElem = setdiff(adjElems, eGlob);
            if isempty(otherElem), continue; end
            otherElem = otherElem(1);
            if eSet(otherElem), continue; end
            if ~any(ismember(otherElem, cat(1, parts.elemIdx))), continue; end

            % Interface edge — outward normal
            xA = node(vaGlob, :);  xB = node(vbGlob, :);
            edgeVec = xB - xA;  L = norm(edgeVec);
            vOpp = setdiff(elem(eGlob, :), [vaGlob, vbGlob]);
            xC = node(vOpp, :);
            edgeMid = (xA + xB) / 2;
            toC = xC - edgeMid;
            nCand = [-edgeVec(2), edgeVec(1)] / L;
            if dot(nCand, toC) > 0, nCand = -nCand; end

            vaLoc = g2l(vaGlob);  vbLoc = g2l(vbGlob);
            nEdges = nEdges + 1;
            eInfo(nEdges).locEdgeNodes = [vaLoc, vbLoc];
            eInfo(nEdges).elemOut = otherElem;
            eInfo(nEdges).normal = nCand;
            eInfo(nEdges).length = L;
        end
    end
    edgeInfo{s} = eInfo;
end

% ---- Subdomain systems ----------------------------------------------------
A_loc = cell(nSub, 1);  b_loc = cell(nSub, 1);
u_loc = cell(nSub, 1);  freeLoc = cell(nSub, 1);  bdLoc = cell(nSub, 1);

for s = 1:nSub
    locNode = parts(s).localNode;
    locElem = parts(s).localElem;
    nLoc = size(locNode, 1);
    A_loc{s} = assembleStiffness2D(locNode, locElem);
    M_s = assembleMass2D(locNode, locElem);
    if isnumeric(f)
        b_loc{s} = M_s * (f * ones(nLoc, 1));
    else
        b_loc{s} = M_s * f(locNode(:,1), locNode(:,2));
    end
    isBd = false(nLoc, 1);
    isBd(ismember(parts(s).nodeIdx, bdNodes)) = true;
    bdLoc{s} = find(isBd);  freeLoc{s} = find(~isBd);
    u_loc{s} = zeros(nLoc, 1);
    if ~isempty(bdLoc{s})
        if isnumeric(uD), u_loc{s}(bdLoc{s}) = uD;
        else, u_loc{s}(bdLoc{s}) = uD(locNode(bdLoc{s},1), locNode(bdLoc{s},2)); end
    end
end

% ---- Iteration -------------------------------------------------------------
convHist = zeros(maxIter, 1);

for iter = 1:maxIter
    u_old = u_loc;

    % Step 1: Overlapping OSM subdomain solves
    for s = 1:nSub
        A_mod = A_loc{s};  b_mod = b_loc{s};
        nLoc = size(A_mod, 1);
        eInfo = edgeInfo{s};

        for ei = 1:length(eInfo)
            vaLoc = eInfo(ei).locEdgeNodes(1);
            vbLoc = eInfo(ei).locEdgeNodes(2);
            L = eInfo(ei).length;  n_i = eInfo(ei).normal;
            T_out = eInfo(ei).elemOut;

            sOutCandidates = find(cellfun(@(e) ismember(T_out, e), {parts.elemIdx}));
            sOutCandidates = setdiff(sOutCandidates, s);
            if isempty(sOutCandidates)
                sOutCandidates = find(cellfun(@(e) ismember(T_out, e), {parts.elemIdx}));
            end
            if isempty(sOutCandidates), continue; end
            sOut = sOutCandidates(1);

            vGlob = elem(T_out, :);
            g2lOut = parts(sOut).global2local;
            uVals = zeros(3, 1);
            for vi = 1:3
                locIdx = g2lOut(vGlob(vi));
                uVals(vi) = u_old{sOut}(locIdx);
            end

            x1 = node(vGlob(1), :);  x2 = node(vGlob(2), :);  x3 = node(vGlob(3), :);
            area2 = (x2(1)-x1(1))*(x3(2)-x1(2)) - (x3(1)-x1(1))*(x2(2)-x1(2));
            gx = [(x2(2)-x3(2))/area2, (x3(2)-x1(2))/area2, (x1(2)-x2(2))/area2];
            gy = [(x3(1)-x2(1))/area2, (x1(1)-x3(1))/area2, (x2(1)-x1(1))/area2];
            gradU = [sum(uVals .* gx'), sum(uVals .* gy')];
            dudn = dot(gradU, n_i);

            vaGlob = parts(s).nodeIdx(vaLoc);
            vbGlob = parts(s).nodeIdx(vbLoc);
            u_va = u_old{sOut}(g2lOut(vaGlob));
            u_vb = u_old{sOut}(g2lOut(vbGlob));
            g_a = dudn + alpha * u_va;
            g_b = dudn + alpha * u_vb;

            A_mod(vaLoc, vaLoc) = A_mod(vaLoc, vaLoc) + alpha * L/3;
            A_mod(vbLoc, vbLoc) = A_mod(vbLoc, vbLoc) + alpha * L/3;
            A_mod(vaLoc, vbLoc) = A_mod(vaLoc, vbLoc) + alpha * L/6;
            A_mod(vbLoc, vaLoc) = A_mod(vbLoc, vaLoc) + alpha * L/6;
            b_mod(vaLoc) = b_mod(vaLoc) + L/3 * g_a + L/6 * g_b;
            b_mod(vbLoc) = b_mod(vbLoc) + L/6 * g_a + L/3 * g_b;
        end

        fL = freeLoc{s};  bL = bdLoc{s};
        if ~isempty(bL)
            u_solve = zeros(nLoc, 1);
            if isnumeric(uD), u_solve(bL) = uD;
            else
                locNd = parts(s).localNode;
                u_solve(bL) = uD(locNd(bL,1), locNd(bL,2));
            end
            b_mod(fL) = b_mod(fL) - A_mod(fL, bL) * u_solve(bL);
            u_solve(fL) = A_mod(fL, fL) \ b_mod(fL);
            u_loc{s} = u_solve;
        else
            u_loc{s} = A_mod \ b_mod;
        end
    end

    % Step 2: Coarse correction
    u_glob = zeros(N, 1);  count = zeros(N, 1);
    for s = 1:nSub
        gIdx = parts(s).nodeIdx;
        u_glob(gIdx) = u_glob(gIdx) + u_loc{s};
        count(gIdx) = count(gIdx) + 1;
    end
    u_glob(count > 0) = u_glob(count > 0) ./ count(count > 0);
    u_free = u_glob(freeNodes);
    r_free = b_f - A_ff * u_free;
    r_H = R_H * r_free;
    if strcmp(coarseSolver, 'chol')
        e_H = L_H \ (L_H' \ r_H);
    else
        e_H = U_H \ (L_H \ (PLU * r_H));
    end
    e_free = P_H * e_H;
    u_free = u_free + e_free;
    u_glob = zeros(N, 1);
    u_glob(freeNodes) = u_free;
    for s = 1:nSub
        gIdx = parts(s).nodeIdx;
        u_loc{s} = u_glob(gIdx);
    end

    % Convergence
    changeNorm = 0;
    for s = 1:nSub
        changeNorm = max(changeNorm, norm(u_loc{s} - u_old{s}, inf));
    end
    convHist(iter) = changeNorm;
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
