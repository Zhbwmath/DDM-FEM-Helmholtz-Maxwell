addpath(genpath(pwd));

[nd, el, bd] = squaremesh([0,1,0,1], 0.25);
[~, edgeIdx, edgeSign] = edgeMesh2D(el);
NE = max(edgeIdx(:));
fprintf('Mesh: N=%d, NT=%d, NE=%d\n', size(nd,1), size(el,1), NE);

A = assembleCurlCurl2D(nd, el);
M = assembleNedMass2D(nd, el);
fprintf('A: %dx%d, nnz=%d, issym=%d\n', size(A,1), size(A,2), nnz(A), issymmetric(A));
fprintf('M: %dx%d, nnz=%d, issym=%d\n', size(M,1), size(M,2), nnz(M), issymmetric(M));

% Check nullspace of A: curl(curl) has gradient fields in nullspace
% For a simply connected domain with n×u=0 on boundary, the kernel 
% should be gradients of H01 functions. Let's check the size.

% Check RHS computation
f_rhs = @(x,y) (pi^2+1)*sin(pi*y);
b = assembleNedRHS2D_debug(nd, el, f_rhs);
fprintf('||b|| = %.4e\n', norm(b));

% Check boundary edges
K = A + M;
bdEdges = findBoundaryEdges2D(el, bd);
freeEdges = setdiff(1:NE, bdEdges)';
fprintf('Total edges: %d, boundary: %d, free: %d\n', NE, length(bdEdges), length(freeEdges));

% Solve
u_f = K(freeEdges, freeEdges) \ b(freeEdges);
uh = zeros(NE, 1); uh(freeEdges) = u_f;
fprintf('||uh|| = %.4e\n', norm(uh));

% Check a specific DOF value against exact
u_ex = @(x,y) [sin(pi*y), zeros(size(x))];
% For an edge at y=0.5 (horizontal edge), the DOF should be ∫ u·t ds ≈ sin(π*0.5)*L = L
% Let's find a horizontal interior edge
for e = 1:min(10,NE)
    verts = edgeMesh2D(el);
    v1 = verts(e,1); v2 = verts(e,2);
    fprintf('Edge %d: v%d(%.2f,%.2f)-v%d(%.2f,%.2f) DOF=%.6f\n',...
        e, v1, nd(v1,1), nd(v1,2), v2, nd(v2,1), nd(v2,2), uh(e));
end
fprintf('(first 10 edges)\n');
exit(0);

function b = assembleNedRHS2D_debug(node, elem, f_rhs)
[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem,1);
[lambda_q, weight] = quadtriangle(2);
nQuad = length(weight);
x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
s1=edgeSign(:,1); s2=edgeSign(:,2); s3=edgeSign(:,3);
b=zeros(NE,1);
for q=1:nQuad
    l=lambda_q(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx = f_rhs(px,py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    c1=weight(q)*area.*(fx.*p1x);
    c2=weight(q)*area.*(fx.*p2x);
    c3=weight(q)*area.*(fx.*p3x);
    b=b+accumarray(edgeIdx(:,1),s1.*c1,[NE,1]);
    b=b+accumarray(edgeIdx(:,2),s2.*c2,[NE,1]);
    b=b+accumarray(edgeIdx(:,3),s3.*c3,[NE,1]);
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
