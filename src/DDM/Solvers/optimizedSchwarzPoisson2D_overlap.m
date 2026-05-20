function [u, convHist] = optimizedSchwarzPoisson2D_overlap(node, elem, bdFlag, f, uD, ...
    parts, alpha, tol, maxIter)
% OPTIMIZEDSCHWARZPOISSON2D_OVERLAP  Overlapping OSM for 2D Poisson.
%
%   Subdomains overlap (δ > 0). Robin condition on ∂Ω_i \ ∂Ω:
%     ∂u_i/∂n_i + α u_i = ∂u_j/∂n_i + α u_j
%
%   Flux uses edge-based one-side gradient from the outward element T_out
%   (the element of Ω_j adjacent to the interface edge, outside Ω_i).
%
%   parts must have overlap > 0.

if nargin < 8, alpha = []; end
if nargin < 9, tol = 1e-6; end
if nargin < 10, maxIter = 100; end

nSub = length(parts);
N = size(node, 1);
NT = size(elem, 1);

if isempty(alpha)
    H = (max(node(:,1)) - min(node(:,1))) / sqrt(nSub);
    alpha = 0.5 * pi / H;
end

bdNodes = getBoundaryNodes2D(elem, bdFlag);

% ---- Precompute outward-facing interface edges for each subdomain ----------
% For each subdomain Ω_i, find edges on ∂Ω_i \ ∂Ω and their outward elements.
% An edge is on ∂Ω_i \ ∂Ω if one adjacent element is in Ω_i, the other is outside
% Ω_i AND inside another subdomain (not outside the global domain).

edgeVP = [2 3; 3 1; 1 2];
% Build edge-to-elements map
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

% For each subdomain, find interface edges and outward elements
edgeInfo = cell(nSub, 1);  % edgeInfo{s} = struct with fields:
                           % .locEdgeNodes: [nEdge x 2] local node indices
                           % .elemIn: element index inside Ω_i (local)
                           % .elemOut: element index of T_out (global)
                           % .normal: [nEdge x 2] outward normal n_i
                           % .length: edge length

for s = 1:nSub
    g2l = parts(s).global2local;
    locElem = parts(s).localElem;
    eIdx = parts(s).elemIdx;
    eSet = false(NT, 1);  eSet(eIdx) = true;

    eInfo = struct('locEdgeNodes', {}, 'elemIn', {}, 'elemOut', {}, ...
                   'normal', {}, 'length', {});
    nEdges = 0;

    for ei = 1:size(locElem, 1)
        eGlob = eIdx(ei);
        for k = 1:3
            vaGlob = elem(eGlob, edgeVP(k, 1));
            vbGlob = elem(eGlob, edgeVP(k, 2));
            key = sprintf('%d,%d', min(vaGlob, vbGlob), max(vaGlob, vbGlob));
            adjElems = edge2elem(key);

            % Find element on the other side
            otherElem = setdiff(adjElems, eGlob);
            if isempty(otherElem), continue; end  % boundary edge
            otherElem = otherElem(1);

            % Is the other element outside Ω_i but in some Ω_j?
            if eSet(otherElem), continue; end  % both sides in Ω_i → interior

            % This edge is on ∂Ω_i \ ∂Ω (or global boundary)
            % Check if it's on global boundary
            if ~any(ismember(otherElem, cat(1, parts.elemIdx)))
                continue;  % global boundary, not a DDM interface
            end

            % Found interface edge. Compute outward normal n_i
            % n_i points FROM element eGlob (in Ω_i) toward otherElem
            % The edge goes from va to vb. The normal is the outward normal
            % of element eGlob across this edge.
            % For triangle with edge opposite vertex v_opp:
            % outward normal points from v_opp toward the edge midpoint

            % Simpler: compute edge tangent and outward normal by geometry
            xA = node(vaGlob, :);  xB = node(vbGlob, :);
            edgeVec = xB - xA;  L = norm(edgeVec);
            % Outward normal from eGlob: cross edge tangent
            % Two possible normals: ±(edgeVec_perp) / L
            % Pick the one pointing away from the third vertex of eGlob

            % Third vertex of eGlob (not on this edge)
            vOpp = setdiff(elem(eGlob, :), [vaGlob, vbGlob]);
            xC = node(vOpp, :);
            edgeMid = (xA + xB) / 2;
            toC = xC - edgeMid;

            % Candidate normal (rotate edgeVec 90° CCW)
            nCand = [-edgeVec(2), edgeVec(1)] / L;

            % If nCand points toward xC (the third vertex), flip it
            if dot(nCand, toC) > 0
                nCand = -nCand;
            end
            % Now nCand points OUTWARD from eGlob

            % Local node indices
            vaLoc = g2l(vaGlob);  vbLoc = g2l(vbGlob);

            nEdges = nEdges + 1;
            eInfo(nEdges).locEdgeNodes = [vaLoc, vbLoc];
            eInfo(nEdges).elemIn = ei;  % local element index in Ω_i
            eInfo(nEdges).elemOut = otherElem;  % global element index
            eInfo(nEdges).normal = nCand;
            eInfo(nEdges).length = L;
        end
    end
    edgeInfo{s} = eInfo;
    fprintf('  Sub %d: %d interface edges on dOmega\\partialOmega\n', s, nEdges);
end

% ---- Pre-build subdomain systems -------------------------------------------
A_loc = cell(nSub, 1);  M_loc = cell(nSub, 1);
b_loc = cell(nSub, 1);  u_loc = cell(nSub, 1);
freeLoc = cell(nSub, 1);  bdLoc = cell(nSub, 1);
locNodeCell = cell(nSub, 1);

for s = 1:nSub
    locNode = parts(s).localNode;
    locElem = parts(s).localElem;
    nLoc = size(locNode, 1);
    locNodeCell{s} = locNode;

    A_loc{s} = assembleStiffness2D(locNode, locElem);
    M_loc{s} = assembleMass2D(locNode, locElem);

    if isnumeric(f)
        b_loc{s} = M_loc{s} * (f * ones(nLoc, 1));
    else
        b_loc{s} = M_loc{s} * f(locNode(:,1), locNode(:,2));
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
end

% ---- Iteration -------------------------------------------------------------
convHist = zeros(maxIter, 1);

for iter = 1:maxIter
    u_old = u_loc;

    % Solve each subdomain
    for s = 1:nSub
        A_mod = A_loc{s};
        b_mod = b_loc{s};
        nLoc = size(A_mod, 1);

        % Add Robin terms on ∂Ω_i \ ∂Ω
        eInfo = edgeInfo{s};
        for ei = 1:length(eInfo)
            vaLoc = eInfo(ei).locEdgeNodes(1);
            vbLoc = eInfo(ei).locEdgeNodes(2);
            L = eInfo(ei).length;
            n_i = eInfo(ei).normal;
            T_out = eInfo(ei).elemOut;  % global element index

            % Compute ∇u on T_out (P1: piecewise constant)
            % Find which subdomain T_out belongs to (prefer neighbor ≠ s)
            sOutCandidates = find(cellfun(@(e) ismember(T_out, e), {parts.elemIdx}));
            sOutCandidates = setdiff(sOutCandidates, s);
            if isempty(sOutCandidates)
                sOutCandidates = find(cellfun(@(e) ismember(T_out, e), {parts.elemIdx}));
            end
            if isempty(sOutCandidates), continue; end
            sOut = sOutCandidates(1);

            % Global vertices of T_out
            vGlob = elem(T_out, :);
            % Get u values
            g2lOut = parts(sOut).global2local;
            uVals = zeros(3, 1);
            for vi = 1:3
                locIdx = g2lOut(vGlob(vi));
                uVals(vi) = u_old{sOut}(locIdx);
            end

            % Element geometry for T_out (P1 gradient)
            x1 = node(vGlob(1), :);  x2 = node(vGlob(2), :);  x3 = node(vGlob(3), :);
            area2 = (x2(1)-x1(1))*(x3(2)-x1(2)) - (x3(1)-x1(1))*(x2(2)-x1(2));

            % ∇φ_i on T_out
            gx = [(x2(2)-x3(2))/area2, (x3(2)-x1(2))/area2, (x1(2)-x2(2))/area2];
            gy = [(x3(1)-x2(1))/area2, (x1(1)-x3(1))/area2, (x2(1)-x1(1))/area2];

            % ∇u on T_out
            gradU = [sum(uVals .* gx'), sum(uVals .* gy')];

            % du_j/dn_i = ∇u_j · n_i
            dudn = dot(gradU, n_i);

            % Robin data at edge endpoints
            u_a = u_old{sOut}(g2lOut(vGlob(1)));  % use actual vertex values
            u_b = u_old{sOut}(g2lOut(vGlob(2)));
            % Actually need values at edge endpoints, not element vertices
            % For the two edge endpoint vertices:
            vaGlob = parts(s).nodeIdx(vaLoc);
            vbGlob = parts(s).nodeIdx(vbLoc);
            u_va = u_old{sOut}(g2lOut(vaGlob));
            u_vb = u_old{sOut}(g2lOut(vbGlob));

            g_a = dudn + alpha * u_va;
            g_b = dudn + alpha * u_vb;

            % Edge mass: (L/6) * [2 1; 1 2]
            % Add to stiffness: α * M_e
            A_mod(vaLoc, vaLoc) = A_mod(vaLoc, vaLoc) + alpha * L/3;
            A_mod(vbLoc, vbLoc) = A_mod(vbLoc, vbLoc) + alpha * L/3;
            A_mod(vaLoc, vbLoc) = A_mod(vaLoc, vbLoc) + alpha * L/6;
            A_mod(vbLoc, vaLoc) = A_mod(vbLoc, vaLoc) + alpha * L/6;

            % Add to RHS: M_e * g
            b_mod(vaLoc) = b_mod(vaLoc) + L/3 * g_a + L/6 * g_b;
            b_mod(vbLoc) = b_mod(vbLoc) + L/6 * g_a + L/3 * g_b;
        end

        % Dirichlet BC and solve
        fL = freeLoc{s};  bL = bdLoc{s};
        if ~isempty(bL)
            u_solve = zeros(nLoc, 1);
            if isnumeric(uD)
                u_solve(bL) = uD;
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

    % Convergence: max change
    changeNorm = 0;
    for s = 1:nSub
        changeNorm = max(changeNorm, norm(u_loc{s} - u_old{s}, inf));
    end
    convHist(iter) = changeNorm;
    if convHist(iter) < tol
        convHist = convHist(1:iter);  break;
    end
end

% Global assembly (average at shared nodes)
u = zeros(N, 1);  count = zeros(N, 1);
for s = 1:nSub
    gIdx = parts(s).nodeIdx;
    u(gIdx) = u(gIdx) + u_loc{s};
    count(gIdx) = count(gIdx) + 1;
end
mask = count > 0;  u(mask) = u(mask) ./ count(mask);
end
