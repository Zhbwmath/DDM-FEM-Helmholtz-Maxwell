addpath(genpath(pwd));
[nd, el, bd] = squaremesh([0,1,0,1], 0.25);
[~, edgeIdx, edgeSign] = edgeMesh2D(el);
NE = max(edgeIdx(:));
NT = size(el,1);

A = assembleCurlCurl2D(nd, el);
M = assembleNedMass2D(nd, el);

% Compute interpolant DOFs correctly via line integrals
u_dof = zeros(NE,1);
u_dof_count = zeros(NE,1);

for t = 1:NT
    v = nd(el(t,:),:);
    for k = 1:3
        switch k
            case 1, va=2; vb=3;
            case 2, va=3; vb=1;
            case 3, va=1; vb=2;
        end
        eid = edgeIdx(t,k);
        sig = edgeSign(t,k);
        a = v(va,:); b = v(vb,:);
        L = norm(b-a);
        tvec = (b-a)/L;
        % Line integral ∫ u·t ds using 3-point Gauss
        [xi,w1d] = gauss1D_mapped(3);
        val = 0;
        for qi = 1:3
            s = (xi(qi)+1)/2;
            pt = a + s*(b-a);
            uval = sin(pi*pt(2));  % x-component of u_exact; y-comp is 0
            val = val + w1d(qi)/2 * L * uval * tvec(1);
        end
        u_dof(eid) = u_dof(eid) + sig * val;
        u_dof_count(eid) = u_dof_count(eid) + 1;
    end
end
% Each interior edge appears twice; average...
u_dof = u_dof ./ max(u_dof_count,1);

fprintf('First 5 exact DOFs: %.6f %.6f %.6f %.6f %.6f\n', u_dof(1:5));

% Now check ||A*u_dof + M*u_dof - b||
f_rhs = @(x,y) (pi^2+1)*sin(pi*y);
[lambda_q, weight] = quadtriangle(4); nQuad = length(weight);
x1=nd(el(:,1),:); x2=nd(el(:,2),:); x3=nd(el(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
s1=edgeSign(:,1); s2=edgeSign(:,2); s3=edgeSign(:,3);

b = zeros(NE,1);
for q=1:nQuad
    l=lambda_q(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py); fy=zeros(size(fx));
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    c1=weight(q)*area.*(fx.*p1x+fy.*p1y);
    c2=weight(q)*area.*(fx.*p2x+fy.*p2y);
    c3=weight(q)*area.*(fx.*p3x+fy.*p3y);
    b=b+accumarray(edgeIdx(:,1),s1.*c1,[NE,1]);
    b=b+accumarray(edgeIdx(:,2),s2.*c2,[NE,1]);
    b=b+accumarray(edgeIdx(:,3),s3.*c3,[NE,1]);
end

r = A*u_dof + M*u_dof - b;
fprintf('||A*u||=%.4e ||M*u||=%.4e ||b||=%.4e ||r||=%.4e\n',...
    norm(A*u_dof), norm(M*u_dof), norm(b), norm(r));
fprintf('ratio = %.2e\n', norm(r)/norm(b));

% Also check: does u_dof satisfy A*u + M*u = b for the FE solution?
bdEdges = findBoundaryEdges2D(el,bd);
freeEdges = setdiff(1:NE,bdEdges)';
K = A + M;
uh_free = K(freeEdges,freeEdges) \ b(freeEdges);
uh = zeros(NE,1); uh(freeEdges)=uh_free;
fprintf('||uh||=%.4e ||u_dof||=%.4e ||uh-u_dof||=%.4e\n',...
    norm(uh), norm(u_dof), norm(uh-u_dof));

% Check: for horizontal interior edge at y=0.5, what's the DOF?
fprintf('\nDOF comparison (first 10):\n');
for e=1:min(10,NE)
    fprintf('  e%d: exact=%.6f  computed=%.6f\n', e, u_dof(e), uh(e));
end

exit(0);

function [x,w]=gauss1D_mapped(n)
switch n
    case 3, x=[-sqrt(3/5);0;sqrt(3/5)]; w=[5/9;8/9;5/9];
end
end

function bdEdges=findBoundaryEdges2D(elem,bdFlag)
[~,edgeIdx]=edgeMesh2D(elem); bdEdges=[];
for k=1:3
    isBd=bdFlag(:,k)==1;
    if any(isBd), bdEdges=[bdEdges;edgeIdx(isBd,k)]; end
end
bdEdges=unique(bdEdges);
end
