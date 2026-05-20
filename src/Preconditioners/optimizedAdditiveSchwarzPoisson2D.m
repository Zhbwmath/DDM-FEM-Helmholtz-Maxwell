function applyPrecon = optimizedAdditiveSchwarzPoisson2D(node, elem, bdFlag, parts, freeNodes, alpha)
% OPTIMIZEDADDITIVESCHWARZPOISSON2D  One-level OAS preconditioner for Poisson.
%
%   Local problems use Dirichlet conditions on the physical boundary and
%   Robin conditions on artificial subdomain boundaries:
%       A_j = K_j + alpha * M_{Gamma_j}.
%
%   The additive correction is weighted by the nodal partition of unity over
%   overlapping subdomains.

if nargin < 6 || isempty(alpha)
    alpha = pi;
end

Nf = length(freeNodes);
global2reduced = zeros(max(freeNodes), 1);
global2reduced(freeNodes) = (1:Nf)';

nSub = length(parts);
nodeCount = zeros(size(node, 1), 1);
for s = 1:nSub
    nodeCount(parts(s).nodeIdx) = nodeCount(parts(s).nodeIdx) + 1;
end

physBdNodes = getBoundaryNodes2D(elem, bdFlag);
physBdEdges = physicalBoundaryEdges2D(elem, bdFlag);

locSolvers = cell(nSub, 1);
locRedIdx = cell(nSub, 1);
locWeights = cell(nSub, 1);

for s = 1:nSub
    gIdx = parts(s).nodeIdx;
    locFree = find(~ismember(gIdx, physBdNodes));
    freeGlobal = gIdx(locFree);
    freeGlobal = intersect(freeGlobal, freeNodes, 'stable');

    if isempty(freeGlobal)
        continue;
    end

    redPos = global2reduced(freeGlobal);
    localPos = parts(s).global2local(freeGlobal);

    Kloc = assembleStiffness2D(parts(s).localNode, parts(s).localElem);
    artFlag = artificialBoundaryFlags2D(parts(s), physBdEdges);
    Mbloc = assembleBoundaryMass2D(parts(s).localNode, parts(s).localElem, artFlag);
    Aloc = Kloc + alpha * Mbloc;
    Aloc = Aloc(localPos, localPos);

    try
        locSolvers{s} = chol(Aloc);
    catch
        [L, U, P] = lu(Aloc);
        locSolvers{s} = {L, U, P};
    end
    locRedIdx{s} = redPos;
    locWeights{s} = 1 ./ nodeCount(freeGlobal);
end

    function x = applyImpl(r)
        x = zeros(Nf, 1);
        for j = 1:nSub
            redPos = locRedIdx{j};
            if isempty(redPos), continue; end

            solver = locSolvers{j};
            rj = r(redPos);
            if ismatrix(solver) && size(solver, 1) == size(solver, 2)
                zj = solver \ (solver' \ rj);
            else
                zj = solver{2} \ (solver{1} \ (solver{3} * rj));
            end
            x(redPos) = x(redPos) + locWeights{j} .* zj;
        end
    end

applyPrecon = @applyImpl;
end


function bdEdges = physicalBoundaryEdges2D(elem, bdFlag)
edgeVerts = [2 3; 3 1; 1 2];
bdEdges = zeros(0, 2);
for k = 1:3
    isBd = bdFlag(:, k) == 1;
    if any(isBd)
        bdEdges = [bdEdges; sort(elem(isBd, edgeVerts(k, :)), 2)]; %#ok<AGROW>
    end
end
bdEdges = unique(bdEdges, 'rows');
end


function artFlag = artificialBoundaryFlags2D(part, physBdEdges)
locElem = part.localElem;
NT = size(locElem, 1);
edgeVerts = [2 3; 3 1; 1 2];

allEdges = zeros(3*NT, 2);
for k = 1:3
    rows = (k-1)*NT + (1:NT);
    allEdges(rows, :) = sort(locElem(:, edgeVerts(k, :)), 2);
end
[~, ~, ic] = unique(allEdges, 'rows');
counts = accumarray(ic, 1);

artFlag = zeros(NT, 3);
for k = 1:3
    rows = (k-1)*NT + (1:NT);
    isSubBoundary = counts(ic(rows)) == 1;

    globEdges = sort(part.nodeIdx(locElem(:, edgeVerts(k, :))), 2);
    isPhysical = ismember(globEdges, physBdEdges, 'rows');
    artFlag(:, k) = isSubBoundary & ~isPhysical;
end
end
