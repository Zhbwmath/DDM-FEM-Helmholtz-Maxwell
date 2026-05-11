addpath(genpath(pwd));

% Single element: reference triangle (0,0), (1,0), (0,1)
node = [0,0; 1,0; 0,1];
elem = [1,2,3];
[~, edgeIdx, edgeSign] = edgeMesh2D(elem);

A = assembleCurlCurl2D(node, elem);
M = assembleNedMass2D(node, elem);

fprintf('Single element - NE_1 local matrices:\n');
fprintf('Stiffness (curl-curl):\n'); disp(full(A));
fprintf('Mass:\n'); disp(full(M));

% Expected values for reference triangle:
% |T| = 0.5
% Curls: curl1=2(∇λ₂×∇λ₃)=2, curl2=2(∇λ₃×∇λ₁)=2, curl3=2(∇λ₁×∇λ₂)=2
% K_loc(i,j) = |T|*curl_i*curl_j = 0.5*2*2 = 2 for all i,j
fprintf('Expected stiffness: all entries = 2\n');
fprintf('Expected mass (diag) ≈ 1/6, (off-diag) ≈ -1/8 or so\n');

% Now manually compute mass for one entry
x1=node(1,:); x2=node(2,:); x3=node(3,:);
area2=(x2(1)-x1(1))*(x3(2)-x1(2))-(x3(1)-x1(1))*(x2(2)-x1(2));
area=abs(area2)/2;
g1=[(x2(2)-x3(2))/area2,(x3(1)-x2(1))/area2];
g2=[(x3(2)-x1(2))/area2,(x1(1)-x3(1))/area2];
g3=[(x1(2)-x2(2))/area2,(x2(1)-x1(1))/area2];
fprintf('\nGradients:\n g1=[%.4f,%.4f]\n g2=[%.4f,%.4f]\n g3=[%.4f,%.4f]\n',...
    g1(1),g1(2),g2(1),g2(2),g3(1),g3(2));

% Curl (scalar) = 2*(g_i × g_j)
curl1=2*(g2(1)*g3(2)-g2(2)*g3(1));
curl2=2*(g3(1)*g1(2)-g3(2)*g1(1));
curl3=2*(g1(1)*g2(2)-g1(2)*g2(1));
fprintf('Curls: c1=%.4f c2=%.4f c3=%.4f\n', curl1, curl2, curl3);
fprintf('|T|*c1^2 = %.4f (expect 2)\n', area*curl1^2);

% Mass with 3-pt quadrature
[lambda_q,w]=quadtriangle(2);
M_manual=zeros(3,3);
for q=1:length(w)
    l=lambda_q(q,:);
    p1=[l(2)*g3(1)-l(3)*g2(1), l(2)*g3(2)-l(3)*g2(2)];
    p2=[l(3)*g1(1)-l(1)*g3(1), l(3)*g1(2)-l(1)*g3(2)];
    p3=[l(1)*g2(1)-l(2)*g1(1), l(1)*g2(2)-l(2)*g1(2)];
    M_manual(1,1)=M_manual(1,1)+w(q)*area*dot(p1,p1);
    M_manual(1,2)=M_manual(1,2)+w(q)*area*dot(p1,p2);
    M_manual(1,3)=M_manual(1,3)+w(q)*area*dot(p1,p3);
    M_manual(2,2)=M_manual(2,2)+w(q)*area*dot(p2,p2);
    M_manual(2,3)=M_manual(2,3)+w(q)*area*dot(p2,p3);
    M_manual(3,3)=M_manual(3,3)+w(q)*area*dot(p3,p3);
end
M_manual(2,1)=M_manual(1,2); M_manual(3,1)=M_manual(1,3); M_manual(3,2)=M_manual(2,3);
fprintf('\nManual mass:\n'); disp(M_manual);
fprintf('Code mass:\n'); disp(full(M));
exit(0);
