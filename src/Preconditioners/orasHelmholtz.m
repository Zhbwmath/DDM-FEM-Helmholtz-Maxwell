function applyPrecon = orasHelmholtz(node, elem, bdFlag, k, parts, degree, solverMode, useParfor)
% ORASHELMHOLTZ  ORAS preconditioner for Helmholtz.
%   B_h^{-1} = Σ_j R̃_{h,j} A_{h,j}^{-1} R_{h,j}
%
%   applyPrecon = ORASHELMHOLTZ(node, elem, bdFlag, k, parts)
%   applyPrecon = ORASHELMHOLTZ(node, elem, bdFlag, k, parts, degree)
%   applyPrecon = ORASHELMHOLTZ(..., degree, solverMode)
%   applyPrecon = ORASHELMHOLTZ(..., degree, solverMode, useParfor)
%
%   solverMode: 'lu' (default) — sparse LU with row/column permutations
%               'direct'     — A \ b each time, slower per iteration, less memory
%   useParfor:  false (default) — sequential for loop
%               true — uses parfor (requires active parallel pool)
%
%   Partition of unity: χ_j(node) = 1 / (#subdomains containing node)
%   A_{h,j} = K_j - k² M_j - ik Mb_j  (impedance on ALL of ∂Ω_j)

if nargin < 6, degree = 1; end
if nargin < 7, solverMode = 'lu'; end
if nargin < 8, useParfor = false; end

if degree == 1
    nodeH = node;
    elemH = elem;
else
    [nodeH, elemH] = extendMesh2D(node, elem, degree);
end

N = size(nodeH,1);  nSub = length(parts);
edgeVP = [2 3; 3 1; 1 2];

% ---- Partition of unity weights -------------------------------------------
useWeightFun = isfield(parts, 'weightFun') && ~isempty(parts(1).weightFun);
nodeWeight = zeros(N,1);
for j = 1:nSub
    idxH = unique(elemH(parts(j).elemIdx, :));
    if useWeightFun
        raw = parts(j).weightFun(nodeH(idxH,1), nodeH(idxH,2));
        nodeWeight(idxH) = nodeWeight(idxH) + max(raw(:), 0);
    else
        nodeWeight(idxH) = nodeWeight(idxH) + 1;
    end
end
nodeWeight(nodeWeight == 0) = 1;

% ---- Subdomain matrices + solvers -----------------------------------------
locSolvers = cell(nSub,1);
gIdx      = cell(nSub,1);
wgt       = cell(nSub,1);

if useParfor
    parfor j = 1:nSub
        [locSolvers{j}, gIdx{j}, wgt{j}] = setupSubdomain(j, parts, elem, ...
            elemH, nodeH, edgeVP, degree, k, solverMode, nodeWeight, useWeightFun);
    end
else
    for j = 1:nSub
        [locSolvers{j}, gIdx{j}, wgt{j}] = setupSubdomain(j, parts, elem, ...
            elemH, nodeH, edgeVP, degree, k, solverMode, nodeWeight, useWeightFun);
    end
end

% ---- Preconditioner application -------------------------------------------
    function x = applyImpl(r)
        x = zeros(N,1);
        for j = 1:nSub
            rj = r(gIdx{j});
            S  = locSolvers{j};
            if strcmpi(solverMode, 'direct')
                zj = S \ rj;                      % backslash (UMFPACK)
            else
                zj = zeros(size(rj));
                zj(S{4}) = S{2} \ (S{1} \ rj(S{3}));
            end
            x(gIdx{j}) = x(gIdx{j}) + zj .* wgt{j};
        end
    end
applyPrecon = @applyImpl;
end

% ---- Subdomain setup (extracted for parfor compatibility) ------------------
function [solver, gIdx_j, wgt_j] = setupSubdomain(j, parts, elem, ...
    elemH, nodeH, edgeVP, degree, k, solverMode, nodeWeight, useWeightFun)

nE   = length(parts(j).elemIdx);
eIdx = parts(j).elemIdx;

gIdx_j = unique(elemH(eIdx, :));
g2l = zeros(size(nodeH, 1), 1);
g2l(gIdx_j) = (1:length(gIdx_j))';
ln = nodeH(gIdx_j, :);
le = g2l(elemH(eIdx, :));

K = assembleStiffness2D(ln, le, degree);
M = assembleMass2D(ln, le, degree);

% Vectorized boundary edge detection
subEdges  = [elem(eIdx, [edgeVP(1,1), edgeVP(1,2)])
             elem(eIdx, [edgeVP(2,1), edgeVP(2,2)])
             elem(eIdx, [edgeVP(3,1), edgeVP(3,2)])];
subEdgesS = sort(subEdges, 2);
[~, ~, subEid] = unique(subEdgesS, 'rows');
eCount = accumarray(subEid, 1);
bdPerEdge = eCount(subEid) == 1;
bdFlagLoc = reshape(bdPerEdge, [nE, 3]);

Mb    = assembleBoundaryMass2D(ln, le, bdFlagLoc, degree);
A_loc = K - k^2*M - 1i*k*Mb;

if strcmpi(solverMode, 'direct')
    solver = A_loc;
else
    [L, U, p, q] = lu(A_loc, 'vector');
    solver = {L, U, p(:), q(:)};
end

if useWeightFun
    raw = max(parts(j).weightFun(nodeH(gIdx_j,1), nodeH(gIdx_j,2)), 0);
    wgt_j = raw(:) ./ nodeWeight(gIdx_j);
else
    wgt_j = 1 ./ nodeWeight(gIdx_j);
end
end
