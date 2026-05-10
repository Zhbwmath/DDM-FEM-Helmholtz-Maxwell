function Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree)
% ASSEMBLEBOUNDARYMASS2D  Assemble the Pk boundary mass matrix on a 2D mesh.
%
%   Mb_ij = \int_{\partial\Omega} \phi_i \phi_j  ds
%
%   Mb = ASSEMBLEBOUNDARYMASS2D(node, elem, bdFlag)        % default: P1
%   Mb = ASSEMBLEBOUNDARYMASS2D(node, elem, bdFlag, degree) % P1, P2, or P3
%
%   bdFlag(t,k) = 1 if edge k of element t is on the Dirichlet boundary.
%
%   P1 uses closed-form  L/6 * [2 1; 1 2]  on each boundary edge.
%   P2/P3 use 1D Gauss-Legendre quadrature on boundary edges.

if nargin < 4, degree = 1; end

if degree == 1
    Mb = assembleBoundaryMass2D_P1(node, elem, bdFlag);
else
    Mb = assembleBoundaryMass2D_quad(node, elem, bdFlag, degree);
end
end


function Mb = assembleBoundaryMass2D_P1(node, elem, bdFlag)
N = size(node, 1);
maxBd = 3 * size(elem, 1);
ii = zeros(maxBd * 4, 1);
jj = zeros(maxBd * 4, 1);
ss = zeros(maxBd * 4, 1);
idx = 0;

edgeVertex = [2 3; 3 1; 1 2];

for k = 1:3
    bdEdges = (bdFlag(:,k) == 1);
    if ~any(bdEdges), continue; end

    e = elem(bdEdges, :);
    vA = e(:, edgeVertex(k,1));
    vB = e(:, edgeVertex(k,2));
    L = sqrt((node(vB,1)-node(vA,1)).^2 + (node(vB,2)-node(vA,2)).^2);

    nBd = length(L);
    nxt = idx + 1;  idx = idx + 4*nBd;

    ii(nxt:4:idx)   = vA;  jj(nxt:4:idx)   = vA;  ss(nxt:4:idx)   = L/3;
    ii(nxt+1:4:idx) = vB;  jj(nxt+1:4:idx) = vB;  ss(nxt+1:4:idx) = L/3;
    ii(nxt+2:4:idx) = vA;  jj(nxt+2:4:idx) = vB;  ss(nxt+2:4:idx) = L/6;
    ii(nxt+3:4:idx) = vB;  jj(nxt+3:4:idx) = vA;  ss(nxt+3:4:idx) = L/6;
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function Mb = assembleBoundaryMass2D_quad(node, elem, bdFlag, degree)
% Quadrature-based boundary mass for P2/P3 on 2D boundary edges.
% The restriction of a Pk 2D basis to an edge is the 1D Pk basis.

% Extend mesh if needed
if size(elem, 2) == 3
    [node, elem, edge, edgeSign] = extendMesh2D(node, elem, degree);
else
    [~, ~, edge, edgeSign] = extendMesh2D(node, elem(:,1:3), degree);
end

N = size(node, 1);
nLB = size(elem, 2);                     % 6 for P2, 10 for P3

% 1D Gauss quadrature on [0,1]
quadOrder = 2 * degree;
[xi, w1d] = gauss1D_01(quadOrder);
nQuad = length(w1d);

% 1D Lagrange basis on [0,1] at Gauss points
phi1d = lagrange1D(degree, xi);          % nQuad x (degree+1)

% Edge-to-element DOF mapping per local edge
% P2: edge k has DOFs: [vertex_a, vertex_b, edge_midpoint]
% P3: edge k has DOFs: [vertex_a, vertex_b, edge_1/3, edge_2/3]

edgeDofs2D = getEdgeDofs2D(degree);

maxBd = 3 * size(elem, 1);
nEntries = maxBd * (degree+1)^2 * 2;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

edgeVerts = [2 3; 3 1; 1 2];             % edges opposite each vertex

for k = 1:3                               % loop over 3 local edges
    bdEdges = (bdFlag(:,k) == 1);
    if ~any(bdEdges), continue; end

    eBd = elem(bdEdges, :);
    nBd = sum(bdEdges);
    eBdIdx = find(bdEdges);

    % Edge endpoints
    vA = eBd(:, edgeVerts(k,1));
    vB = eBd(:, edgeVerts(k,2));
    L = sqrt((node(vB,1)-node(vA,1)).^2 + (node(vB,2)-node(vA,2)).^2);

    % DOF indices on this edge
    dofIdx = edgeDofs2D(k, :, degree);

    for q = 1:nQuad
        phi_q = phi1d(q, :)';            % (deg+1) x 1
        for a = 1:(degree+1)
            ia = dofIdx(a);
            for b = a:(degree+1)
                ib = dofIdx(b);
                s = w1d(q) * L * (phi_q(a) * phi_q(b));
                nxt = idx + 1;  idx = idx + nBd;
                ii(nxt:idx) = eBd(:, ia);
                jj(nxt:idx) = eBd(:, ib);
                ss(nxt:idx) = s;
                if a ~= b
                    nxt2 = idx + 1;  idx = idx + nBd;
                    ii(nxt2:idx) = eBd(:, ib);
                    jj(nxt2:idx) = eBd(:, a);
                    ss(nxt2:idx) = s;
                end
            end
        end
    end
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


% ===========================================================================
function dofIdx = getEdgeDofs2D(edgeK, degree)
% Return the column indices in elem_e for DOFs on local edge edgeK.
%   edgeK = 1 (opp v1): edge (v2,v3)
%   edgeK = 2 (opp v2): edge (v3,v1)
%   edgeK = 3 (opp v3): edge (v1,v2)
%
%   P2 ordering: v1,v2,v3, e12,e23,e31
%     edge1 (v2,v3): DOFs [2, 3, 5]   (v2, v3, e23)
%     edge2 (v3,v1): DOFs [3, 1, 6]   (v3, v1, e31)
%     edge3 (v1,v2): DOFs [1, 2, 4]   (v1, v2, e12)
%
%   P3 ordering: v1,v2,v3, e12a,e12b, e23a,e23b, e31a,e31b, centroid
%     edge1 (v2,v3): DOFs [2, 3, 6, 7]      (v2, v3, e23a, e23b)
%     edge2 (v3,v1): DOFs [3, 1, 8, 9]      (v3, v1, e31a, e31b)
%     edge3 (v1,v2): DOFs [1, 2, 4, 5]      (v1, v2, e12a, e12b)

switch degree
    case 2
        switch edgeK
            case 1, dofIdx = [2, 3, 5];    % edge (v2,v3)
            case 2, dofIdx = [3, 1, 6];    % edge (v3,v1)
            case 3, dofIdx = [1, 2, 4];    % edge (v1,v2)
        end
    case 3
        switch edgeK
            case 1, dofIdx = [2, 3, 6, 7]; % edge (v2,v3): near v2 then v3
            case 2, dofIdx = [3, 1, 8, 9]; % edge (v3,v1): near v3 then v1
            case 3, dofIdx = [1, 2, 4, 5]; % edge (v1,v2): near v1 then v2
        end
    otherwise
        error('getEdgeDofs2D: degree %d not supported', degree);
end
end


% ===========================================================================
function [x, w] = gauss1D_01(n)
% Gauss-Legendre quadrature on [0,1].
[x_ref, w_ref] = gauss1D_mapped(n);
x = (x_ref + 1) / 2;
w = w_ref / 2;
end


function [x, w] = gauss1D_mapped(n)
% Gauss-Legendre on [-1, 1] (same as in quadtet.m helper).
switch n
    case 1,  x = 0;  w = 2;
    case 2
        x = [-sqrt(1/3); sqrt(1/3)];
        w = [1; 1];
    case 3
        x = [-sqrt(3/5); 0; sqrt(3/5)];
        w = [5/9; 8/9; 5/9];
    case 4
        x = [-sqrt(3/7 + 2/7*sqrt(6/5));
             -sqrt(3/7 - 2/7*sqrt(6/5));
              sqrt(3/7 - 2/7*sqrt(6/5));
              sqrt(3/7 + 2/7*sqrt(6/5))];
        w = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];
    case 5
        x = [-sqrt(5+2*sqrt(10/7))/3;
             -sqrt(5-2*sqrt(10/7))/3;
              0;
              sqrt(5-2*sqrt(10/7))/3;
              sqrt(5+2*sqrt(10/7))/3];
        w = [(322-13*sqrt(70))/900;
             (322+13*sqrt(70))/900;
              128/225;
             (322+13*sqrt(70))/900;
             (322-13*sqrt(70))/900];
    case 6
        x = [-0.932469514203152; -0.661209386466265; -0.238619186083197;
              0.238619186083197;  0.661209386466265;  0.932469514203152];
        w = [0.171324492379170; 0.360761573048139; 0.467913934572691;
             0.467913934572691; 0.360761573048139; 0.171324492379170];
    otherwise
        error('gauss1D: n=%d not implemented', n);
end
end


% ===========================================================================
function phi = lagrange1D(degree, x)
% 1D Lagrange basis on [0,1] at points x (nQuad x 1).
% Returns phi: nQuad x (degree+1).
% Nodes: P1: 0,1  P2: 0,0.5,1  P3: 0,1/3,2/3,1

nQuad = length(x);
switch degree
    case 1
        phi = [(1-x), x];
    case 2
        phi = [(2*x-1).*(x-1), 4*x.*(1-x), x.*(2*x-1)];
    case 3
        phi = zeros(nQuad, 4);
        phi(:,1) = -4.5 * (x - 1/3) .* (x - 2/3) .* (x - 1);
        phi(:,2) = 13.5 * (x - 0) .* (x - 2/3) .* (x - 1);
        phi(:,3) = -13.5 * (x - 0) .* (x - 1/3) .* (x - 1);
        phi(:,4) = 4.5 * (x - 0) .* (x - 1/3) .* (x - 2/3);
    otherwise
        error('lagrange1D: degree %d not supported', degree);
end
end
