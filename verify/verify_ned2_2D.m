% VERIFY_NED2_2D  Convergence test for NE_2 in 2D.
%   u = [y(1-y), 0] on [0,1]^2, n×u=0 on boundary.

fprintf('========== 2D NE_2 Convergence Study ==========\n\n');

u_exact = @(x,y) deal(y.*(1-y), zeros(size(x)));
curl_exact = @(x,y) 2*y - 1;
f_rhs = @(x,y) 2 + y.*(1-y);

nRefine = 3;
fmt = '%-8s  %-8s  %-12s  %-8s  %-12s  %-8s\n';
fprintf(fmt, 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_Hcurl', 'rateHc');
fprintf('%s\n', repmat('-',1,66));

for k = 1:nRefine
    hk = 2^(-k-1);
    [nd, el, bd] = squaremesh([0,1,0,1], hk);

    [~, eidx, es] = edgeMesh2D(el);
    NE = max(eidx(:)); NT = size(el,1);
    Ntot = 2*NE + 2*NT;

    A = assembleNed2CurlCurl2D(nd, el);
    M = assembleNed2Mass2D(nd, el);
    K = A + M;

    % RHS
    b = assembleNed2RHS2D(nd, el, f_rhs);

    % Boundary edges → their 2 DOFs are constrained
    bdE = findBoundaryEdges2D(el, bd);
    bdDOFs = sort([2*(bdE-1)+1; 2*(bdE-1)+2]);
    freeDOFs = setdiff(1:Ntot, bdDOFs)';

    uh = zeros(Ntot,1);
    uh(freeDOFs) = K(freeDOFs, freeDOFs) \ b(freeDOFs);

    [eL2, eHc] = computeNed2Error2D(nd, el, uh, u_exact, curl_exact);

    if k > 1
        rL2 = log(eL2/eL2p)/log(hk/hp);
        rHc = log(eHc/eHcp)/log(hk/hp);
        fprintf(fmt, sprintf('%.4f',hk), sprintf('%d',Ntot), ...
            sprintf('%.4e',eL2), sprintf('%.2f',rL2), sprintf('%.4e',eHc), sprintf('%.2f',rHc));
    else
        fprintf(fmt, sprintf('%.4f',hk), sprintf('%d',Ntot), ...
            sprintf('%.4e',eL2), '-', sprintf('%.4e',eHc), '-');
    end
    eL2p=eL2; eHcp=eHc; hp=hk;
end
fprintf('\nExpected for this NE_2 variant: L2~O(h), H(curl)~O(h)\n');
fprintf('========== Done ==========\n');


function b = assembleNed2RHS2D(node, elem, f_rhs)
[~, eidx, es] = edgeMesh2D(elem);
NE = max(eidx(:)); NT = size(elem,1); Ntot = 2*NE+2*NT; nLocal = 8;
[lambda_q, weight] = quadtriangle(4);
x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
localEdgeCols=[2,3,1]; localEdgeBases=[3,5,1];
gIdx=zeros(NT,nLocal); gSign=zeros(NT,nLocal);
for kk=1:3
    col=localEdgeCols(kk); eidk=eidx(:,col); sigk=es(:,col); b0=localEdgeBases(kk);
    gIdx(:,b0)=2*(eidk-1)+1; gIdx(:,b0+1)=2*(eidk-1)+2;
    gSign(:,b0)=sigk; gSign(:,b0+1)=1;
end
for t=1:NT, gIdx(t,7)=2*NE+2*(t-1)+1; gIdx(t,8)=2*NE+2*(t-1)+2; end; gSign(:,7:8)=1;
b=zeros(Ntot,1);
for q=1:length(weight)
    l=lambda_q(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py); fy=zeros(size(fx));
    phix=zeros(NT,nLocal); phiy=zeros(NT,nLocal);
    phix(:,1)=l(1)*g2(:,1)-l(2)*g1(:,1); phiy(:,1)=l(1)*g2(:,2)-l(2)*g1(:,2);
    c12=l(1)-l(2); phix(:,2)=c12.*phix(:,1); phiy(:,2)=c12.*phiy(:,1);
    phix(:,3)=l(2)*g3(:,1)-l(3)*g2(:,1); phiy(:,3)=l(2)*g3(:,2)-l(3)*g2(:,2);
    c23=l(2)-l(3); phix(:,4)=c23.*phix(:,3); phiy(:,4)=c23.*phiy(:,3);
    phix(:,5)=l(3)*g1(:,1)-l(1)*g3(:,1); phiy(:,5)=l(3)*g1(:,2)-l(1)*g3(:,2);
    c31=l(3)-l(1); phix(:,6)=c31.*phix(:,5); phiy(:,6)=c31.*phiy(:,5);
    phix(:,7)=l(1).*l(2).*g3(:,1); phiy(:,7)=l(1).*l(2).*g3(:,2);
    phix(:,8)=l(2).*l(3).*g1(:,1); phiy(:,8)=l(2).*l(3).*g1(:,2);
    for kk=1:nLocal
        c=2*weight(q)*area.*(fx.*phix(:,kk)+fy.*phiy(:,kk));
        b=b+accumarray(gIdx(:,kk),gSign(:,kk).*c,[Ntot,1]);
    end
end
end

function bdE = findBoundaryEdges2D(elem, bdFlag)
[~, eidx] = edgeMesh2D(elem);
bdFlag_to_e = [2,3,1]; bdE=[];
for k=1:3, isBd=bdFlag(:,k)==1; if any(isBd), bdE=[bdE;eidx(isBd,bdFlag_to_e(k))]; end; end
bdE=unique(bdE);
end

function [eL2, eHc] = computeNed2Error2D(node, elem, uh, u_exact, curl_exact)
[~, eidx, es] = edgeMesh2D(elem);
NE=max(eidx(:)); NT=size(elem,1); nLocal=8;
[lambda_q, weight]=quadtriangle(6);
x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
localEdgeCols=[2,3,1]; localEdgeBases=[3,5,1];
gIdx=zeros(NT,nLocal); gSign=zeros(NT,nLocal);
for kk=1:3
    col=localEdgeCols(kk); eidk=eidx(:,col); sigk=es(:,col); b0=localEdgeBases(kk);
    gIdx(:,b0)=2*(eidk-1)+1; gIdx(:,b0+1)=2*(eidk-1)+2;
    gSign(:,b0)=sigk; gSign(:,b0+1)=1;
end
for t=1:NT, gIdx(t,7)=2*NE+2*(t-1)+1; gIdx(t,8)=2*NE+2*(t-1)+2; end; gSign(:,7:8)=1;
eL2=0; eHc=0;
for q=1:length(weight)
    l=lambda_q(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    [uex,uey]=u_exact(px,py); curlex=curl_exact(px,py);
    phix=zeros(NT,nLocal); phiy=zeros(NT,nLocal); curll=zeros(NT,nLocal);
    phix(:,1)=l(1)*g2(:,1)-l(2)*g1(:,1); phiy(:,1)=l(1)*g2(:,2)-l(2)*g1(:,2);
    c12=l(1)-l(2); phix(:,2)=c12.*phix(:,1); phiy(:,2)=c12.*phiy(:,1);
    phix(:,3)=l(2)*g3(:,1)-l(3)*g2(:,1); phiy(:,3)=l(2)*g3(:,2)-l(3)*g2(:,2);
    c23=l(2)-l(3); phix(:,4)=c23.*phix(:,3); phiy(:,4)=c23.*phiy(:,3);
    phix(:,5)=l(3)*g1(:,1)-l(1)*g3(:,1); phiy(:,5)=l(3)*g1(:,2)-l(1)*g3(:,2);
    c31=l(3)-l(1); phix(:,6)=c31.*phix(:,5); phiy(:,6)=c31.*phiy(:,5);
    phix(:,7)=l(1).*l(2).*g3(:,1); phiy(:,7)=l(1).*l(2).*g3(:,2);
    phix(:,8)=l(2).*l(3).*g1(:,1); phiy(:,8)=l(2).*l(3).*g1(:,2);
    curll(:,1)=2*(g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1));
    curll(:,2)=2*c12.*curll(:,1)+2*phix(:,1).*(g1(:,2)-g2(:,2))-2*phiy(:,1).*(g1(:,1)-g2(:,1));
    curll(:,3)=2*(g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1));
    curll(:,4)=2*c23.*curll(:,3)+2*phix(:,3).*(g2(:,2)-g3(:,2))-2*phiy(:,3).*(g2(:,1)-g3(:,1));
    curll(:,5)=2*(g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1));
    curll(:,6)=2*c31.*curll(:,5)+2*phix(:,5).*(g3(:,2)-g1(:,2))-2*phiy(:,5).*(g3(:,1)-g1(:,1));
    curll(:,7)=(l(2)*g1(:,1)+l(1)*g2(:,1)).*g3(:,2)-(l(2)*g1(:,2)+l(1)*g2(:,2)).*g3(:,1);
    curll(:,8)=(l(3)*g2(:,1)+l(2)*g3(:,1)).*g1(:,2)-(l(3)*g2(:,2)+l(2)*g3(:,2)).*g1(:,1);
    uhx=zeros(NT,1); uhy=zeros(NT,1); cuh=zeros(NT,1);
    for kk=1:nLocal
        uv=uh(gIdx(:,kk)); ss=gSign(:,kk);
        uhx=uhx+ss.*uv.*phix(:,kk); uhy=uhy+ss.*uv.*phiy(:,kk);
        cuh=cuh+ss.*uv.*curll(:,kk);
    end
    ex=uhx-uex; ey=uhy-uey; ec=cuh-curlex;
    wa=2*weight(q)*area;
    eL2=eL2+sum(wa.*(ex.^2+ey.^2));
    eHc=eHc+sum(wa.*(ex.^2+ey.^2+ec.^2));
end
eL2=sqrt(eL2); eHc=sqrt(eHc);
end
