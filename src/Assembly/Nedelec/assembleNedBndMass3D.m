function Mb = assembleNedBndMass3D(node, elem, bdFlag)
% ASSEMBLENEDBNDMASS3D  Tangential boundary mass for NE_1 in 3D.
%
%   Mb_ij = \int_{\partial\Omega} (n×φ_i) · (n×φ_j)  dS
%
%   Uses 2D quadrature on boundary faces.  On each boundary face, only
%   the 3 edges of that face contribute non-zero tangential traces.

[~, edgeIdx, edgeSign] = edgeMesh3D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

% Face definitions: face k is opposite vertex k
faceVerts = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
% Edges on each face (local edge indices)
faceEdges = {[4,6,5], [2,6,3], [1,5,3], [1,4,2]};
% Local edge pairs: (1,2)=k1, (1,3)=k2, (1,4)=k3, (2,3)=k4, (2,4)=k5, (3,4)=k6

[lambda2d, w2d] = quadtriangle(2);       % 2D quadrature on faces
nQuad = length(w2d);

% For each boundary face, use 2D quadrature
bdVals_i = []; bdVals_j = []; bdVals_s = [];

for f = 1:4
    isFace = (bdFlag(:,f) == 1);
    if ~any(isFace), continue; end
    idx = find(isFace);
    nBd = length(idx);

    fv = faceVerts{f};
    fe = faceEdges{f};                    % 3 local edge indices for this face

    % Face vertices
    a = node(elem(idx, fv(1)), :);
    b = node(elem(idx, fv(2)), :);
    c = node(elem(idx, fv(3)), :);

    % Face area and normal (outward)
    ab = b - a;  ac = c - a;
    cr = cross(ab, ac);
    area = 0.5 * sqrt(cr(:,1).^2 + cr(:,2).^2 + cr(:,3).^2);
    n = cr ./ (2 * area);                 % unit normal (outward)

    % Barycentric gradients of the tet
    v1n=node(elem(idx,1),:); v2n=node(elem(idx,2),:);
    v3n=node(elem(idx,3),:); v4n=node(elem(idx,4),:);
    e12t=v2n-v1n; e13t=v3n-v1n; e14t=v4n-v1n;
    detJt=e12t(:,1).*(e13t(:,2).*e14t(:,3)-e13t(:,3).*e14t(:,2)) ...
         +e12t(:,2).*(e13t(:,3).*e14t(:,1)-e13t(:,1).*e14t(:,3)) ...
         +e12t(:,3).*(e13t(:,1).*e14t(:,2)-e13t(:,2).*e14t(:,1));
    invJt=1./detJt;
    g2t=cross(e13t,e14t).*invJt; g3t=cross(e14t,e12t).*invJt; g4t=cross(e12t,e13t).*invJt;
    g1t=-(g2t+g3t+g4t);
    Gt=cell(4,1); Gt{1}=g1t; Gt{2}=g2t; Gt{3}=g3t; Gt{4}=g4t;

    edges_all = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

    for q = 1:nQuad
        l = lambda2d(q,:);
        % Barycentric coordinates on the tet for this face point
        % Face f is opposite vertex f, so λ_f = 0 on this face
        lam = zeros(nBd, 4);
        for k_local = 1:3
            lam(:, fv(k_local)) = l(k_local);
        end

        % Evaluate basis φ_ij at this point
        phi = zeros(nBd, 6, 3);
        for k = 1:6
            i = edges_all(k,1); j = edges_all(k,2);
            gi = Gt{i}; gj = Gt{j};
            phi(:,k,1) = lam(:,i).*gj(:,1) - lam(:,j).*gi(:,1);
            phi(:,k,2) = lam(:,i).*gj(:,2) - lam(:,j).*gi(:,2);
            phi(:,k,3) = lam(:,i).*gj(:,3) - lam(:,j).*gi(:,3);
        end

        % Tangential trace: n × φ (vector)
        % (n×φ)·(n×ψ) = (φ·ψ) - (φ·n)(ψ·n)  [for unit n]
        % Or compute directly: n×φ
        for ei = 1:3                    % 3 edges of this face
            ke = fe(ei);                 % local edge index (1..6)
            for ej = 1:3
                le = fe(ej);
                % n×φ_ke · n×φ_le
                % = (n·n)(φ_ke·φ_le) - (φ_ke·n)(φ_le·n) = (φ_ke·φ_le) - (φ_ke·n)(φ_le·n)
                phi_dot = phi(:,ke,1).*phi(:,le,1) + phi(:,ke,2).*phi(:,le,2) + phi(:,ke,3).*phi(:,le,3);
                phi_ndot_ke = phi(:,ke,1).*n(:,1) + phi(:,ke,2).*n(:,2) + phi(:,ke,3).*n(:,3);
                phi_ndot_le = phi(:,le,1).*n(:,1) + phi(:,le,2).*n(:,2) + phi(:,le,3).*n(:,3);
                val = phi_dot - phi_ndot_ke .* phi_ndot_le;

                % Global DOF indices with signs
                eid_ke = edgeIdx(idx, ke);  sig_ke = edgeSign(idx, ke);
                eid_le = edgeIdx(idx, le);  sig_le = edgeSign(idx, le);

                s = 2 * w2d(q) * area .* sig_ke .* sig_le .* val;
                bdVals_i = [bdVals_i; eid_ke];
                bdVals_j = [bdVals_j; eid_le];
                bdVals_s = [bdVals_s; s];
            end
        end
    end
end

Mb = sparse(bdVals_i, bdVals_j, bdVals_s, NE, NE);
end
