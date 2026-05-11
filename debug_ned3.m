addpath(genpath(pwd));
[nd, el, bd] = squaremesh([0,1,0,1], 0.25);
[~, edgeIdx, edgeSign] = edgeMesh2D(el);
NE = max(edgeIdx(:));

A = assembleCurlCurl2D(nd, el);
M = assembleNedMass2D(nd, el);
K = A + M;

% Compute exact DOFs: u_e = \int_edge u_exact · t ds
u_exact = @(x,y) [sin(pi*y), zeros(size(x))];
u_dof_exact = zeros(NE, 1);

for t = 1:size(el,1)
    v = nd(el(t,:),:);
    % Edge 1: v2-v3
    L = norm(v(3,:)-v(2,:)); mid = (v(2,:)+v(3,:))/2;
    u_dof_exact(edgeIdx(t,1)) = edgeSign(t,1) * L * dot(u_exact(mid(1),mid(2)), (v(3,:)-v(2,:))/L);
    % Edge 2: v3-v1
    L = norm(v(1,:)-v(3,:)); mid = (v(3,:)+v(1,:))/2;
    u_dof_exact(edgeIdx(t,2)) = edgeSign(t,2) * L * dot(u_exact(mid(1),mid(2)), (v(1,:)-v(3,:))/L);
    % Edge 3: v1-v2
    L = norm(v(2,:)-v(1,:)); mid = (v(1,:)+v(2,:))/2;
    u_dof_exact(edgeIdx(t,3)) = edgeSign(t,3) * L * dot(u_exact(mid(1),mid(2)), (v(2,:)-v(1,:))/L);
end

% Check: does K*u_dof_exact ≈ RHS?
f_rhs = @(x,y) (pi^2+1)*sin(pi*y);
% Simple RHS: b_i = \int f·φ_i. For u_exact, f = (π²+1)u, and 
% \int f·φ_i = (π²+1) \int u·φ_i
% Since u·φ_i is related to DOFs...
% Actually, check property: for NE_1 edge basis φ_e,
% ∫ u_h · φ_e dx = Σ_f M_ef u_f (where M is mass matrix)
% So K*u = (A+M)*u should equal RHS vector b = ∫ f·φ

% Let me compute b properly
[lambda_q, weight] = quadtriangle(2); nQuad = length(weight);
x1=nd(el(:,1),:); x2=nd(el(:,2),:); x3=nd(el(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
s1=edgeSign(:,1); s2=edgeSign(:,2); s3=edgeSign(:,3);

b = zeros(NE,1);
for q = 1:nQuad
    l = lambda_q(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx = f_rhs(px,py); fy = zeros(size(fx));
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

% Check residual
r = K * u_dof_exact - b;
fprintf('||K*u_ex - b|| = %.4e, ||b|| = %.4e, ratio = %.4e\n', norm(r), norm(b), norm(r)/norm(b));

% Check DOF vs exact solution
% For horizontal edge at y=0.5
fprintf('\nExact DOFs for first 10 edges:\n');
for e = 1:min(10, NE)
    [v1,v2] = deal(edgeMesh2D(el));  % need the edge list
end
% Let me just check the DOF values
fprintf('Exact first 5 DOFs: %.4f %.4f %.4f %.4f %.4f\n', u_dof_exact(1:5));

exit(0);
