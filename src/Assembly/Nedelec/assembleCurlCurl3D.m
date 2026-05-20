function A = assembleCurlCurl3D(node, elem)
% ASSEMBLECURLCURL3D  Assemble the NE_1 curl-curl stiffness matrix in 3D.
%
%   A_ij = \int_\Omega curl(φ_i) · curl(φ_j)  dV
%
%   NE_1 has constant curl (3D vector) per element → no quadrature needed.

[~, edgeIdx, edgeSign] = edgeMesh3D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);
nLocal = 6;

% Barycentric gradients (NT x 4 x 3)
v1 = node(elem(:,1),:); v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:); v4 = node(elem(:,4),:);
e12=v2-v1; e13=v3-v1; e14=v4-v1;

detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ)/6;
invJ = 1./detJ;

% ∇λ_2 = cross(e13,e14)/detJ, ∇λ_3 = cross(e14,e12)/detJ, ∇λ_4 = cross(e12,e13)/detJ
% ∇λ_1 = -(∇λ_2 + ∇λ_3 + ∇λ_4)
g2 = cross(e13, e14) .* invJ;  % NT x 3
g3 = cross(e14, e12) .* invJ;
g4 = cross(e12, e13) .* invJ;
g1 = -(g2 + g3 + g4);

% Pack into gradLam: gradLam(:,i,:) = ∇λ_i
gradLam = zeros(NT, 4, 3);
gradLam(:,1,:) = g1; gradLam(:,2,:) = g2;
gradLam(:,3,:) = g3; gradLam(:,4,:) = g4;

% Curl of each local basis: curl(φ_ij) = 2 ∇λ_i × ∇λ_j
edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
c = zeros(NT, nLocal, 3);
for k = 1:nLocal
    i = edges(k,1);  j = edges(k,2);
    gi = squeeze(gradLam(:,i,:));  % NT x 3
    gj = squeeze(gradLam(:,j,:));
    c(:,k,1) = 2*(gi(:,2).*gj(:,3) - gi(:,3).*gj(:,2));
    c(:,k,2) = 2*(gi(:,3).*gj(:,1) - gi(:,1).*gj(:,3));
    c(:,k,3) = 2*(gi(:,1).*gj(:,2) - gi(:,2).*gj(:,1));
end

% Local stiffness: K_loc(p,q) = |T| * c_p · c_q
Kloc = zeros(NT, nLocal, nLocal);
for p = 1:nLocal
    for q = p:nLocal
        Kloc(:,p,q) = volume .* (c(:,p,1).*c(:,q,1) + c(:,p,2).*c(:,q,2) + c(:,p,3).*c(:,q,3));
        if p ~= q, Kloc(:,q,p) = Kloc(:,p,q); end
    end
end

% ---- Sparse assembly --------------------------------------------------
nEntries = NT * nLocal * nLocal;
ii = zeros(nEntries,1); jj = zeros(nEntries,1); ss = zeros(nEntries,1);
idx = 0;
for p = 1:nLocal
    ep = edgeIdx(:,p);  sp = edgeSign(:,p);
    for q = 1:nLocal
        eq = edgeIdx(:,q);  sq = edgeSign(:,q);
        s = sp .* sq .* Kloc(:,p,q);
        nxt = idx+1; idx = idx+NT;
        ii(nxt:idx)=ep; jj(nxt:idx)=eq; ss(nxt:idx)=s;
    end
end
A = sparse(ii, jj, ss, NE, NE);
end
