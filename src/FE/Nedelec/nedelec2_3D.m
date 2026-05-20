function [phi, curl_phi] = nedelec2_3D(lambda, gradLambda)
% NEDELEC2_3D  Evaluate NE_2 basis on a 3D tetrahedron (20 DOFs).
%
%   [phi, curl_phi] = NEDELEC2_3D(lambda, gradLambda)
%
%   Local DOF ordering (20):
%     1-12: edge DOFs (2 per edge, 6 edges)
%     13-20: face DOFs (2 per face, 4 faces)

nQuad = size(lambda, 1);
G = cell(4,1);
for i = 1:4, G{i} = gradLambda(i,:); end

phi = zeros(nQuad, 20, 3);
curl_phi = zeros(nQuad, 20, 3);

edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
faceBubbles = {[2,3,4; 3,4,2],[1,4,3; 4,3,1],[1,2,4; 2,4,1],[1,3,2; 3,2,1]};

% ---- Edge basis (DOFs 1-12) -----------------------------------------------
for k = 1:6
    i = edges(k,1); j = edges(k,2);
    gi = G{i}; gj = G{j};
    li = lambda(:,i); lj = lambda(:,j);
    d0 = 2*(k-1) + 1;  d1 = d0 + 1;

    % DOF 0: phi = lam_i*grad(lam_j) - lam_j*grad(lam_i), odd parity
    phix = li*gj(1) - lj*gi(1);
    phiy = li*gj(2) - lj*gi(2);
    phiz = li*gj(3) - lj*gi(3);
    phi(:,d0,1)=phix; phi(:,d0,2)=phiy; phi(:,d0,3)=phiz;
    curl_phi(:,d0,1) = 2*(gi(2)*gj(3) - gi(3)*gj(2));
    curl_phi(:,d0,2) = 2*(gi(3)*gj(1) - gi(1)*gj(3));
    curl_phi(:,d0,3) = 2*(gi(1)*gj(2) - gi(2)*gj(1));

    % DOF 1: phi = (lam_i-lam_j)*(lam_i*grad(lam_j) - lam_j*grad(lam_i)), even parity
    cij = li - lj;
    phi(:,d1,1)=cij.*phix; phi(:,d1,2)=cij.*phiy; phi(:,d1,3)=cij.*phiz;
    dgi = gi - gj;
    curl_phi(:,d1,1) = (dgi(2).*phiz-dgi(3).*phiy) + cij.*curl_phi(:,d0,1);
    curl_phi(:,d1,2) = (dgi(3).*phix-dgi(1).*phiz) + cij.*curl_phi(:,d0,2);
    curl_phi(:,d1,3) = (dgi(1).*phiy-dgi(2).*phix) + cij.*curl_phi(:,d0,3);
end

% ---- Face basis (DOFs 13-20) ----------------------------------------------
for f = 1:4
    fb = faceBubbles{f};
    d0 = 12 + 2*(f-1) + 1;  d1 = d0 + 1;

    % Face basis B_abc = lam_a*lam_b*grad(lam_c) - lam_a*lam_c*grad(lam_b)
    a=fb(1,1); b=fb(1,2); c=fb(1,3);
    ga=G{a}; gb=G{b}; gc=G{c};
    la=lambda(:,a); lb=lambda(:,b); lc=lambda(:,c);
    phi(:,d0,1)=la.*lb*gc(1)-la.*lc*gb(1);
    phi(:,d0,2)=la.*lb*gc(2)-la.*lc*gb(2);
    phi(:,d0,3)=la.*lb*gc(3)-la.*lc*gb(3);
    gx1=lb.*ga(1)+la.*gb(1); gy1=lb.*ga(2)+la.*gb(2); gz1=lb.*ga(3)+la.*gb(3);
    gx2=lc.*ga(1)+la.*gc(1); gy2=lc.*ga(2)+la.*gc(2); gz2=lc.*ga(3)+la.*gc(3);
    curl_phi(:,d0,1)=(gy1*gc(3)-gz1*gc(2))-(gy2*gb(3)-gz2*gb(2));
    curl_phi(:,d0,2)=(gz1*gc(1)-gx1*gc(3))-(gz2*gb(1)-gx2*gb(3));
    curl_phi(:,d0,3)=(gx1*gc(2)-gy1*gc(1))-(gx2*gb(2)-gy2*gb(1));

    % Face basis B_abc = lam_a*lam_b*grad(lam_c) - lam_a*lam_c*grad(lam_b)
    a2=fb(2,1); b2=fb(2,2); c2=fb(2,3);
    ga2=G{a2}; gb2=G{b2}; gc2=G{c2};
    la2=lambda(:,a2); lb2=lambda(:,b2); lc2=lambda(:,c2);
    phi(:,d1,1)=la2.*lb2*gc2(1)-la2.*lc2*gb2(1);
    phi(:,d1,2)=la2.*lb2*gc2(2)-la2.*lc2*gb2(2);
    phi(:,d1,3)=la2.*lb2*gc2(3)-la2.*lc2*gb2(3);
    hx1=lb2.*ga2(1)+la2.*gb2(1); hy1=lb2.*ga2(2)+la2.*gb2(2); hz1=lb2.*ga2(3)+la2.*gb2(3);
    hx2=lc2.*ga2(1)+la2.*gc2(1); hy2=lc2.*ga2(2)+la2.*gc2(2); hz2=lc2.*ga2(3)+la2.*gc2(3);
    curl_phi(:,d1,1)=(hy1*gc2(3)-hz1*gc2(2))-(hy2*gb2(3)-hz2*gb2(2));
    curl_phi(:,d1,2)=(hz1*gc2(1)-hx1*gc2(3))-(hz2*gb2(1)-hx2*gb2(3));
    curl_phi(:,d1,3)=(hx1*gc2(2)-hy1*gc2(1))-(hx2*gb2(2)-hy2*gb2(1));
end
end
