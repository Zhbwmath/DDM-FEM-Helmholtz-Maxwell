function M = assembleNedMass3D(node, elem)
% ASSEMBLENEDMASS3D  Assemble the NE_1 mass matrix in 3D.
%
%   M_ij = \int_\Omega φ_i · φ_j  dV
%
%   Uses quadrature on the reference tetrahedron (order 2).

[~, edgeIdx, edgeSign] = edgeMesh3D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);
nLocal = 6;

[lambda_q, weight] = quadtet(2);
nQuad = length(weight);

v1=node(elem(:,1),:); v2=node(elem(:,2),:);
v3=node(elem(:,3),:); v4=node(elem(:,4),:);
e12=v2-v1; e13=v3-v1; e14=v4-v1;
detJ=e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
    +e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
    +e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume=abs(detJ)/6; invJ=1./detJ;
g2=cross(e13,e14).*invJ; g3=cross(e14,e12).*invJ; g4=cross(e12,e13).*invJ;
g1=-(g2+g3+g4);
G=cell(4,1); G{1}=g1; G{2}=g2; G{3}=g3; G{4}=g4;
edges=[1 2;1 3;1 4;2 3;2 4;3 4];

nEntries=NT*nLocal*nLocal*nQuad;
ii=zeros(nEntries,1); jj=zeros(nEntries,1); ss=zeros(nEntries,1);
idx=0;

for q=1:nQuad
    l=lambda_q(q,:);
    phi_q=zeros(NT,nLocal,3);
    for k=1:nLocal
        i=edges(k,1); j=edges(k,2);
        gi=G{i}; gj=G{j};
        phi_q(:,k,1)=l(i)*gj(:,1)-l(j)*gi(:,1);
        phi_q(:,k,2)=l(i)*gj(:,2)-l(j)*gi(:,2);
        phi_q(:,k,3)=l(i)*gj(:,3)-l(j)*gi(:,3);
    end
    w=6*weight(q)*volume;                % Jacobian: |T|/(1/6) = 6|T|
    for p=1:nLocal
        ep=edgeIdx(:,p); sp=edgeSign(:,p);
        for qq=1:nLocal
            eqq=edgeIdx(:,qq); sqq=edgeSign(:,qq);
            dp=phi_q(:,p,1).*phi_q(:,qq,1)+phi_q(:,p,2).*phi_q(:,qq,2)+phi_q(:,p,3).*phi_q(:,qq,3);
            s=sp.*sqq.*w.*dp;
            nxt=idx+1; idx=idx+NT;
            ii(nxt:idx)=ep; jj(nxt:idx)=eqq; ss(nxt:idx)=s;
        end
    end
end
M=sparse(ii(1:idx),jj(1:idx),ss(1:idx),NE,NE);
end
