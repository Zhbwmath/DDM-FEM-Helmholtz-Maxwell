function Mb = assembleBoundaryMass3D(node, elem, bdFlag, degree)
% ASSEMBLEBOUNDARYMASS3D  Assemble the Pk boundary mass matrix on a 3D mesh.
%
%   Mb_ij = \int_{\partial\Omega} \phi_i \phi_j  ds
%
%   Mb = ASSEMBLEBOUNDARYMASS3D(node, elem, bdFlag)        % default: P1
%   Mb = ASSEMBLEBOUNDARYMASS3D(node, elem, bdFlag, degree) % P1, P2, or P3
%
%   The restriction of a Pk 3D basis to a triangular face is the Pk 2D basis.

if nargin < 4, degree = 1; end

if degree == 1
    Mb = assembleBoundaryMass3D_P1(node, elem, bdFlag);
else
    Mb = assembleBoundaryMass3D_quad(node, elem, bdFlag, degree);
end
end


function Mb = assembleBoundaryMass3D_P1(node, elem, bdFlag)
N = size(node, 1);  NT = size(elem, 1);
faceVerts = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
maxBd = 4*NT;
ii = zeros(maxBd*9,1);  jj = zeros(maxBd*9,1);  ss = zeros(maxBd*9,1);
idx = 0;

for f = 1:4
    bdFaces = (bdFlag(:,f)==1);
    if ~any(bdFaces), continue; end
    fv = faceVerts{f};  e = elem(bdFaces,:);
    vA = e(:,fv(1));  vB = e(:,fv(2));  vC = e(:,fv(3));
    AB = node(vB,:)-node(vA,:);  AC = node(vC,:)-node(vA,:);
    cr = cross(AB, AC);
    area = 0.5 * sqrt(cr(:,1).^2 + cr(:,2).^2 + cr(:,3).^2);

    nBd = length(area);  nxt = idx+1;  idx = idx + 9*nBd;
    ii(nxt:9:idx)=vA; jj(nxt:9:idx)=vA; ss(nxt:9:idx)=area/6;
    ii(nxt+1:9:idx)=vB; jj(nxt+1:9:idx)=vB; ss(nxt+1:9:idx)=area/6;
    ii(nxt+2:9:idx)=vC; jj(nxt+2:9:idx)=vC; ss(nxt+2:9:idx)=area/6;
    ii(nxt+3:9:idx)=vA; jj(nxt+3:9:idx)=vB; ss(nxt+3:9:idx)=area/12;
    ii(nxt+4:9:idx)=vB; jj(nxt+4:9:idx)=vA; ss(nxt+4:9:idx)=area/12;
    ii(nxt+5:9:idx)=vA; jj(nxt+5:9:idx)=vC; ss(nxt+5:9:idx)=area/12;
    ii(nxt+6:9:idx)=vC; jj(nxt+6:9:idx)=vA; ss(nxt+6:9:idx)=area/12;
    ii(nxt+7:9:idx)=vB; jj(nxt+7:9:idx)=vC; ss(nxt+7:9:idx)=area/12;
    ii(nxt+8:9:idx)=vC; jj(nxt+8:9:idx)=vB; ss(nxt+8:9:idx)=area/12;
end
Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function Mb = assembleBoundaryMass3D_quad(node, elem, bdFlag, degree)
% Quadrature-based boundary mass for P2/P3 on 3D tetrahedral faces.
% The restriction of a Pk tetrahedral basis to a face is the Pk 2D basis.

if size(elem, 2) == 4
    [node, elem] = extendMesh3D(node, elem, degree);
end

N = size(node, 1);
nLB = size(elem, 2);

% 2D quadrature on the reference triangle for face integration
quadOrder = 2 * degree;
[lambda2d, w2d] = quadtriangle(quadOrder);
nQuad = length(w2d);

% 2D Lagrange basis on the triangle at these quadrature points
[phi_face, ~] = lagrange2D(degree, lambda2d);
% phi_face: nQuad x nLB2d  where nLB2d = 3, 6, 10 for deg 1,2,3

% DOF indices on each face
faceDofs3D = getFaceDofs3D(degree);
faceVerts = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};

maxBd = 4 * size(elem, 1);
nEntries = maxBd * nQuad * (degree+1)*(degree+2)/2 * 2;
nEntries = min(nEntries, maxBd * 100);   % practical upper bound
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for f = 1:4
    bdFaces = (bdFlag(:,f) == 1);
    if ~any(bdFaces), continue; end

    fv = faceVerts{f};  e = elem(bdFaces,:);
    vA = e(:,fv(1));  vB = e(:,fv(2));  vC = e(:,fv(3));
    AB = node(vB,:)-node(vA,:);  AC = node(vC,:)-node(vA,:);
    cr = cross(AB, AC);
    area = 0.5 * sqrt(cr(:,1).^2 + cr(:,2).^2 + cr(:,3).^2);
    nBd = length(area);

    dofIdx = faceDofs3D(f, :, degree);

    for q = 1:nQuad
        phi_q = phi_face(q, :)';         % nLB2d x 1
        for a_idx = 1:length(dofIdx)
            ia = dofIdx(a_idx);
            for b_idx = a_idx:length(dofIdx)
                ib = dofIdx(b_idx);
                s = w2d(q) * area * (phi_q(a_idx) * phi_q(b_idx));
                nxt = idx + 1;  idx = idx + nBd;
                ii(nxt:idx) = e(:,ia);  jj(nxt:idx) = e(:,ib);  ss(nxt:idx) = s;
                if a_idx ~= b_idx
                    nxt2 = idx + 1;  idx = idx + nBd;
                    ii(nxt2:idx) = e(:,ib);  jj(nxt2:idx) = e(:,ia);  ss(nxt2:idx) = s;
                end
            end
        end
    end
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function dofIdx = getFaceDofs3D(faceK, degree)
% Return column indices in elem_e for DOFs on local face faceK (opp vertex faceK).
%
% P2 (10 nodes): v1..v4,  e12,e13,e14,e23,e24,e34
%   Face 1 (opp v1): v2,v3,v4, e23,e34,e42 → DOFs [2,3,4,8,10,9]
%   Face 2 (opp v2): v1,v3,v4, e13,e34,e41 → DOFs [1,3,4,6,10,7]
%   Face 3 (opp v3): v1,v2,v4, e12,e24,e41 → DOFs [1,2,4,5,9,7]
%   Face 4 (opp v4): v1,v2,v3, e12,e23,e13 → DOFs [1,2,3,5,8,6]
%
% P3 (20 nodes): v1..v4, 6 edges*2 + 4 faces
%   Each face has: 3 vertices + 6 edge points + 1 face centroid = 10 DOFs

switch degree
    case 2
        switch faceK
            case 1, dofIdx = [2,3,4, 8,10,9];  % opp v1
            case 2, dofIdx = [1,3,4, 6,10,7];  % opp v2
            case 3, dofIdx = [1,2,4, 5,9,7];   % opp v3
            case 4, dofIdx = [1,2,3, 5,8,6];   % opp v4
        end
    case 3
        % Face 1 (opp v1): v2,v3,v4 + edges (2,3),(3,4),(4,2) + face centroid f234
        % Edges: (2,3)→11,12, (3,4)→15,16, (4,2)→14,13 (reversed for consistent orientation)
        %   wait, edge (4,2) is edge (2,4) = nodes 13(near v2),14(near v4)
        %   On face opp v1, we need v2→v3→v4 ordering:
        %   edge (2,3): near v2=11, near v3=12
        %   edge (3,4): near v3=15, near v4=16
        %   edge (4,2): near v4=14, near v2=13
        %   face: f234=20
        switch faceK
            case 1  % opp v1: vertices v2,v3,v4
                dofIdx = [2,3,4, 11,12, 15,16, 14,13, 20];
            case 2  % opp v2: vertices v1,v3,v4
                % Edges: (1,3)→7,8  (3,4)→15,16  (4,1)→10,9
                dofIdx = [1,3,4, 7,8, 15,16, 10,9, 19];
            case 3  % opp v3: vertices v1,v2,v4
                % Edges: (1,2)→5,6  (2,4)→13,14  (4,1)→10,9 (reversed: 9,10)
                dofIdx = [1,2,4, 5,6, 13,14, 9,10, 18];
            case 4  % opp v4: vertices v1,v2,v3
                % Edges: (1,2)→5,6  (2,3)→11,12  (3,1)→8,7 (reversed)
                dofIdx = [1,2,3, 5,6, 11,12, 8,7, 17];
        end
    otherwise
        error('getFaceDofs3D: degree %d not supported', degree);
end
end
