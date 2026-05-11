addpath(genpath(pwd));

% Test: on a reference triangle, verify that phi_i · t_j = delta_ij / L_i
% Triangle vertices: (0,0), (1,0), (0,1)
node = [0,0; 1,0; 0,1];
x1=node(1,:); x2=node(2,:); x3=node(3,:);
area2 = (x2(1)-x1(1))*(x3(2)-x1(2)) - (x3(1)-x1(1))*(x2(2)-x1(2));
area = abs(area2)/2;
g1 = [(x2(2)-x3(2))/area2, (x3(1)-x2(1))/area2];
g2 = [(x3(2)-x1(2))/area2, (x1(1)-x3(1))/area2];
g3 = [(x1(2)-x2(2))/area2, (x2(1)-x1(1))/area2];

fprintf('Edge lengths: e1(v2-v3)=%.4f, e2(v3-v1)=%.4f, e3(v1-v2)=%.4f\n',...
    norm(x3-x2), norm(x1-x3), norm(x2-x1));

% Tangents
t1 = (x3-x2)/norm(x3-x2);  % edge 1: v2→v3
t2 = (x1-x3)/norm(x1-x3);  % edge 2: v3→v1
t3 = (x2-x1)/norm(x2-x1);  % edge 3: v1→v2

% Test at several points on each edge
% On edge 3 (v1-v2): lambda = [s, 1-s, 0]
for s = [0, 0.25, 0.5, 0.75, 1.0]
    l = [s, 1-s, 0];
    phi1 = l(2)*g3 - l(3)*g2;  % vector
    phi2 = l(3)*g1 - l(1)*g3;
    phi3 = l(1)*g2 - l(2)*g1;
    fprintf('Edge3 s=%.2f: phi1·t3=%.4f phi2·t3=%.4f phi3·t3=%.4f (expect 0,0,%.4f)\n',...
        s, dot(phi1,t3), dot(phi2,t3), dot(phi3,t3), 1/norm(x2-x1));
end
exit(0);
