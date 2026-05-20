function [u, convHist] = orasSchwarz2D(node, elem, bdFlag, k, f, g, ...
    parts, tol, maxIter, degree)
% ORASSCHWARZ2D  Schwarz iteration with transmission data.
%
%   PDE Schwarz method (Gong et al. 2022, eq. 2.7-2.10):
%     u^{n+1} = Σ_j χ_j u_j^{n+1} with impedance transmission.
%
%   Transmission data: imp(u^n) = ∂u^n/∂n - ik u^n on ∂Ω_j \ ∂Ω,
%   computed from local residual as imp = A_loc*u_loc - b_loc.
if nargin < 9, tol = 1e-6; end
if nargin < 10, maxIter = 200; end
if nargin < 11, degree = 1; end

N = size(node, 1); nSub = length(parts);
edgeVP = [2 3; 3 1; 1 2]; NT = size(elem, 1);

% Edge-to-elements map
edge2elem = containers.Map('KeyType','char','ValueType','any');
for e = 1:NT
    for kv = 1:3
        va = elem(e, edgeVP(kv,1)); vb = elem(e, edgeVP(kv,2));
        key = sprintf('%d,%d', min(va,vb), max(va,vb));
        if isKey(edge2elem, key), edge2elem(key) = [edge2elem(key), e];
        else, edge2elem(key) = e; end
    end
end

% ---- Precompute ----------------------------------------------------------
A_loc = cell(nSub,1); b_loc = cell(nSub,1); M_art = cell(nSub,1);
gIdx = cell(nSub,1); chiV = cell(nSub,1); solver = cell(nSub,1);
nodeCount = zeros(N,1);

for j = 1:nSub
    nodeCount(parts(j).nodeIdx) = nodeCount(parts(j).nodeIdx) + 1;
end

for j = 1:nSub
    ln = parts(j).localNode; le = parts(j).localElem;
    nL = size(ln,1); eIdx = parts(j).elemIdx;

    K = assembleStiffness2D(ln, le, degree);
    M = assembleMass2D(ln, le, degree);

    % Boundary flags: all boundary edges vs artificial only
    bdAll = zeros(length(eIdx),3); bdArt = zeros(length(eIdx),3);
    for ei = 1:length(eIdx)
        eg = eIdx(ei);
        for kv = 1:3
            va = elem(eg, edgeVP(kv,1)); vb = elem(eg, edgeVP(kv,2));
            key = sprintf('%d,%d', min(va,vb), max(va,vb));
            adj = edge2elem(key); other = setdiff(adj, eg);
            if isempty(other) || ~ismember(other(1), eIdx)
                bdAll(ei,kv) = 1;
                if ~isempty(other)
                    for js = 1:nSub
                        if js~=j && ismember(other(1), parts(js).elemIdx)
                            bdArt(ei,kv)=1; break;
                        end
                    end
                end
            end
        end
    end

    MbAll = assembleBoundaryMass2D(ln, le, bdAll, degree);
    MbArt = assembleBoundaryMass2D(ln, le, bdArt, degree);

    A_loc{j} = K - k^2*M - 1i*k*MbAll;
    M_art{j} = MbArt;

    if isnumeric(f), b_loc{j} = M*(f*ones(nL,1));
    else, b_loc{j} = M*f(ln(:,1), ln(:,2)); end
    if nargin>=7 && ~isempty(g)
        if isnumeric(g), b_loc{j}=b_loc{j}+MbAll*(g*ones(nL,1));
        else, b_loc{j}=b_loc{j}+MbAll*g(ln(:,1),ln(:,2)); end
    end

    [L,U,P] = lu(A_loc{j}); solver{j} = {L,U,P};
    gIdx{j} = parts(j).nodeIdx;
    chiV{j} = 1 ./ nodeCount(parts(j).nodeIdx);
end

% ---- Iteration -----------------------------------------------------------
u = zeros(N,1); convHist = zeros(maxIter,1);

for iter = 1:maxIter
    uNewLoc = cell(nSub,1);
    for j = 1:nSub
        uj = u(gIdx{j});  % restrict
        bMod = b_loc{j};

        % Homogeneous impedance on artificial boundary
        % (stationary ORAS preconditioner, eq 2.5)
        % Transmission data = 0 on ∂Ω_j \ ∂Ω

        % Solve
        S = solver{j}; zj = S{2} \ (S{1} \ (S{3} * bMod));
        uNewLoc{j} = zj .* chiV{j};
    end

    uNew = zeros(N,1);
    for j = 1:nSub
        idx = gIdx{j}; val = uNewLoc{j};
        for k = 1:length(idx)
            uNew(idx(k)) = uNew(idx(k)) + val(k);
        end
    end

    chg = norm(uNew - u) / max(norm(uNew),1);
    u = uNew; convHist(iter) = chg;
    if chg < tol, convHist = convHist(1:iter); break; end
end
end
