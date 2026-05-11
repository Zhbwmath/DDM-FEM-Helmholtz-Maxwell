% VERIFY_NED1_3D  Convergence test for NE_1 in 3D.
%   u = [0, 0, x(1-x)y(1-y)]^T on [0,1]^3, n×u=0 on boundary.
%   curl(u) = [x(1-x)(1-2y), -y(1-y)(1-2x), 0]
%   f = curl(curl u) + u

fprintf('========== 3D NE_1 Convergence Study ==========\n\n');

u_exact = @(x,y,z) deal(zeros(size(x)), zeros(size(x)), x.*(1-x).*y.*(1-y));
curl_exact = @(x,y,z) deal(x.*(1-x).*(1-2*y), -y.*(1-y).*(1-2*x), zeros(size(x)));
f_rhs = @(x,y,z) deal(zeros(size(x)), zeros(size(x)), 2*y.*(1-y)+2*x.*(1-x)+x.*(1-x).*y.*(1-y));

nRefine = 3;
fmt = '%-6s  %-8s  %-12s  %-8s  %-12s  %-8s\n';
fprintf(fmt, 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_Hcurl', 'rateHc');
fprintf('%s\n', repmat('-',1,66));

for k = 1:nRefine
    hk = 2^(-k-1);
    [nd, el, bd] = cubemesh([0,1,0,1,0,1], hk);
    [~, eidx_tmp] = edgeMesh3D(el); NE = max(eidx_tmp(:));
    NT = size(el,1);

    A = assembleCurlCurl3D(nd, el);
    M = assembleNedMass3D(nd, el);
    K = A + M;

    % RHS
    [~, eidx, es] = edgeMesh3D(el);
    b = assembleNedRHS3D(nd, el, f_rhs);

    % Boundary
    bdE = findBoundaryEdges3D(el, bd);
    freeE = setdiff(1:NE, bdE)';
    uh = zeros(NE,1); uh(freeE) = K(freeE,freeE) \ b(freeE);

    % Errors
    [eL2, eHc] = computeNedError3D(nd, el, uh, u_exact, curl_exact);

    if k > 1
        rL2 = log(eL2/eL2p)/log(hk/hp);
        rHc = log(eHc/eHcp)/log(hk/hp);
        fprintf(fmt, sprintf('%.4f',hk), sprintf('%d',NE), ...
            sprintf('%.4e',eL2), sprintf('%.2f',rL2), sprintf('%.4e',eHc), sprintf('%.2f',rHc));
    else
        fprintf(fmt, sprintf('%.4f',hk), sprintf('%d',NE), ...
            sprintf('%.4e',eL2), '-', sprintf('%.4e',eHc), '-');
    end
    eL2p=eL2; eHcp=eHc; hp=hk;
end
fprintf('\nExpected: NE_1 L2~O(h), H(curl)~O(h)\n');
fprintf('========== Done ==========\n');


function b = assembleNedRHS3D(node, elem, f_rhs)
[~, eidx, es] = edgeMesh3D(elem);
NE = max(eidx(:)); NT = size(elem,1);
[lambda_q, weight] = quadtet(2);
v1=node(elem(:,1),:); v2=node(elem(:,2),:); v3=node(elem(:,3),:); v4=node(elem(:,4),:);
e12=v2-v1; e13=v3-v1; e14=v4-v1;
detJ=e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
    +e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
    +e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume=abs(detJ)/6; invJ=1./detJ;
g2=cross(e13,e14).*invJ; g3=cross(e14,e12).*invJ; g4=cross(e12,e13).*invJ;
g1=-(g2+g3+g4); G=cell(4,1); G{1}=g1; G{2}=g2; G{3}=g3; G{4}=g4;
edges=[1 2;1 3;1 4;2 3;2 4;3 4]; nLocal=6;
b=zeros(NE,1);
for q=1:length(weight)
    l=lambda_q(q,:);
    px=l(1)*v1(:,1)+l(2)*v2(:,1)+l(3)*v3(:,1)+l(4)*v4(:,1);
    py=l(1)*v1(:,2)+l(2)*v2(:,2)+l(3)*v3(:,2)+l(4)*v4(:,2);
    pz=l(1)*v1(:,3)+l(2)*v2(:,3)+l(3)*v3(:,3)+l(4)*v4(:,3);
    [fx,fy,fz]=f_rhs(px,py,pz);
    for k=1:nLocal
        i=edges(k,1); j=edges(k,2);
        gi=G{i}; gj=G{j};
        phix=l(i)*gj(:,1)-l(j)*gi(:,1);
        phiy=l(i)*gj(:,2)-l(j)*gi(:,2);
        phiz=l(i)*gj(:,3)-l(j)*gi(:,3);
        c=6*weight(q)*volume.*(fx.*phix+fy.*phiy+fz.*phiz);
        b=b+accumarray(eidx(:,k),es(:,k).*c,[NE,1]);
    end
end
end

function bdE = findBoundaryEdges3D(elem, bdFlag)
[~, eidx] = edgeMesh3D(elem);
faceEdges = {[4,6,5], [2,6,3], [1,5,3], [1,4,2]};
bdE = [];
for f = 1:4
    isF = bdFlag(:,f)==1;
    if any(isF)
        for e_local = faceEdges{f}
            bdE = [bdE; eidx(isF, e_local)];
        end
    end
end
bdE = unique(bdE);
end

function [eL2, eHc] = computeNedError3D(node, elem, uh, u_exact, curl_exact)
[~, eidx, es] = edgeMesh3D(elem);
NT = size(elem,1);
[lambda_q, weight] = quadtet(4);
v1=node(elem(:,1),:); v2=node(elem(:,2),:); v3=node(elem(:,3),:); v4=node(elem(:,4),:);
e12=v2-v1; e13=v3-v1; e14=v4-v1;
detJ=e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
    +e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
    +e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume=abs(detJ)/6; invJ=1./detJ;
g2=cross(e13,e14).*invJ; g3=cross(e14,e12).*invJ; g4=cross(e12,e13).*invJ;
g1=-(g2+g3+g4); G=cell(4,1); G{1}=g1; G{2}=g2; G{3}=g3; G{4}=g4;
edges=[1 2;1 3;1 4;2 3;2 4;3 4]; nLocal=6;
% Curl of each basis
c=zeros(NT,nLocal,3);
for k=1:nLocal
    i=edges(k,1); j=edges(k,2);
    gi=G{i}; gj=G{j};
    c(:,k,1)=2*(gi(:,2).*gj(:,3)-gi(:,3).*gj(:,2));
    c(:,k,2)=2*(gi(:,3).*gj(:,1)-gi(:,1).*gj(:,3));
    c(:,k,3)=2*(gi(:,1).*gj(:,2)-gi(:,2).*gj(:,1));
end
eL2=0; eHc=0;
for q=1:length(weight)
    l=lambda_q(q,:);
    px=l(1)*v1(:,1)+l(2)*v2(:,1)+l(3)*v3(:,1)+l(4)*v4(:,1);
    py=l(1)*v1(:,2)+l(2)*v2(:,2)+l(3)*v3(:,2)+l(4)*v4(:,2);
    pz=l(1)*v1(:,3)+l(2)*v2(:,3)+l(3)*v3(:,3)+l(4)*v4(:,3);
    [uex,uey,uez]=u_exact(px,py,pz);
    [cex,cey,cez]=curl_exact(px,py,pz);
    uhx=zeros(NT,1); uhy=zeros(NT,1); uhz=zeros(NT,1);
    cuhx=zeros(NT,1); cuhy=zeros(NT,1); cuhz=zeros(NT,1);
    for k=1:nLocal
        i=edges(k,1); j=edges(k,2);
        gi=G{i}; gj=G{j};
        phix=l(i)*gj(:,1)-l(j)*gi(:,1);
        phiy=l(i)*gj(:,2)-l(j)*gi(:,2);
        phiz=l(i)*gj(:,3)-l(j)*gi(:,3);
        uv=uh(eidx(:,k)); ss=es(:,k);
        uhx=uhx+ss.*uv.*phix; uhy=uhy+ss.*uv.*phiy; uhz=uhz+ss.*uv.*phiz;
        cuhx=cuhx+ss.*uv.*c(:,k,1); cuhy=cuhy+ss.*uv.*c(:,k,2); cuhz=cuhz+ss.*uv.*c(:,k,3);
    end
    ex=uhx-uex; ey=uhy-uey; ez=uhz-uez;
    ecx=cuhx-cex; ecy=cuhy-cey; ecz=cuhz-cez;
    wa=6*weight(q)*volume;
    eL2=eL2+sum(wa.*(ex.^2+ey.^2+ez.^2));
    eHc=eHc+sum(wa.*(ex.^2+ey.^2+ez.^2+ecx.^2+ecy.^2+ecz.^2));
end
eL2=sqrt(eL2); eHc=sqrt(eHc);
end
