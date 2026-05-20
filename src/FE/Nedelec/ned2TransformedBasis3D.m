function [Bx, By, Bz, Cx, Cy, Cz, volume] = ned2TransformedBasis3D(node, elem, lambda, trans)
% NED2TRANSFORMEDBASIS3D  Evaluate oriented 3D NE_2 basis functions.

NT = size(elem, 1);
nLocal = 20;

[G, volume] = tetGradients(node, elem);
[phix, phiy, phiz, curlx, curly, curlz] = localBasis(lambda, G, NT);

Bx = zeros(NT, nLocal); By = zeros(NT, nLocal); Bz = zeros(NT, nLocal);
Cx = zeros(NT, nLocal); Cy = zeros(NT, nLocal); Cz = zeros(NT, nLocal);

% Edge blocks are diagonal signs.
for p = 1:12
    c = trans(:, p, p);
    Bx(:, p) = c .* phix(:, p);
    By(:, p) = c .* phiy(:, p);
    Bz(:, p) = c .* phiz(:, p);
    Cx(:, p) = c .* curlx(:, p);
    Cy(:, p) = c .* curly(:, p);
    Cz(:, p) = c .* curlz(:, p);
end

% Face blocks are 2x2 orientation transforms.
for f = 1:4
    d0 = 12 + 2*(f-1) + 1;
    d1 = d0 + 1;

    t11 = trans(:, d0, d0);
    t21 = trans(:, d1, d0);
    t12 = trans(:, d0, d1);
    t22 = trans(:, d1, d1);

    Bx(:, d0) = t11.*phix(:, d0) + t21.*phix(:, d1);
    By(:, d0) = t11.*phiy(:, d0) + t21.*phiy(:, d1);
    Bz(:, d0) = t11.*phiz(:, d0) + t21.*phiz(:, d1);
    Cx(:, d0) = t11.*curlx(:, d0) + t21.*curlx(:, d1);
    Cy(:, d0) = t11.*curly(:, d0) + t21.*curly(:, d1);
    Cz(:, d0) = t11.*curlz(:, d0) + t21.*curlz(:, d1);

    Bx(:, d1) = t12.*phix(:, d0) + t22.*phix(:, d1);
    By(:, d1) = t12.*phiy(:, d0) + t22.*phiy(:, d1);
    Bz(:, d1) = t12.*phiz(:, d0) + t22.*phiz(:, d1);
    Cx(:, d1) = t12.*curlx(:, d0) + t22.*curlx(:, d1);
    Cy(:, d1) = t12.*curly(:, d0) + t22.*curly(:, d1);
    Cz(:, d1) = t12.*curlz(:, d0) + t22.*curlz(:, d1);
end
end


function [G, volume] = tetGradients(node, elem)
v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

e12 = v2 - v1;
e13 = v3 - v1;
e14 = v4 - v1;

detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));

volume = abs(detJ) / 6;
invJ = 1 ./ detJ;

g2 = cross(e13, e14) .* invJ;
g3 = cross(e14, e12) .* invJ;
g4 = cross(e12, e13) .* invJ;
g1 = -(g2 + g3 + g4);

G = cell(4, 1);
G{1} = g1; G{2} = g2; G{3} = g3; G{4} = g4;
end


function [phix, phiy, phiz, curlx, curly, curlz] = localBasis(lambda, G, NT)
nLocal = 20;
phix = zeros(NT, nLocal); phiy = zeros(NT, nLocal); phiz = zeros(NT, nLocal);
curlx = zeros(NT, nLocal); curly = zeros(NT, nLocal); curlz = zeros(NT, nLocal);

edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for k = 1:6
    i = edges(k, 1);
    j = edges(k, 2);
    gi = G{i};
    gj = G{j};
    li = lambda(i);
    lj = lambda(j);
    d0 = 2*(k-1) + 1;
    d1 = d0 + 1;

    phix(:, d0) = li*gj(:,1) - lj*gi(:,1);
    phiy(:, d0) = li*gj(:,2) - lj*gi(:,2);
    phiz(:, d0) = li*gj(:,3) - lj*gi(:,3);

    curlx(:, d0) = 2*(gi(:,2).*gj(:,3) - gi(:,3).*gj(:,2));
    curly(:, d0) = 2*(gi(:,3).*gj(:,1) - gi(:,1).*gj(:,3));
    curlz(:, d0) = 2*(gi(:,1).*gj(:,2) - gi(:,2).*gj(:,1));

    cij = li - lj;
    dg = gi - gj;
    phix(:, d1) = cij .* phix(:, d0);
    phiy(:, d1) = cij .* phiy(:, d0);
    phiz(:, d1) = cij .* phiz(:, d0);

    curlx(:, d1) = (dg(:,2).*phiz(:,d0) - dg(:,3).*phiy(:,d0)) + cij.*curlx(:,d0);
    curly(:, d1) = (dg(:,3).*phix(:,d0) - dg(:,1).*phiz(:,d0)) + cij.*curly(:,d0);
    curlz(:, d1) = (dg(:,1).*phiy(:,d0) - dg(:,2).*phix(:,d0)) + cij.*curlz(:,d0);
end

faceBubbles = {
    [2,3,4; 3,4,2]
    [1,4,3; 4,3,1]
    [1,2,4; 2,4,1]
    [1,3,2; 3,2,1]
};

for f = 1:4
    fb = faceBubbles{f};
    d0 = 12 + 2*(f-1) + 1;
    d1 = d0 + 1;

    [phix(:, d0), phiy(:, d0), phiz(:, d0), ...
        curlx(:, d0), curly(:, d0), curlz(:, d0)] = faceBasis(fb(1,:), lambda, G);

    [phix(:, d1), phiy(:, d1), phiz(:, d1), ...
        curlx(:, d1), curly(:, d1), curlz(:, d1)] = faceBasis(fb(2,:), lambda, G);
end
end


function [bx, by, bz, cx, cy, cz] = faceBasis(idx, lambda, G)
a = idx(1);
b = idx(2);
c = idx(3);

ga = G{a};
gb = G{b};
gc = G{c};
la = lambda(a);
lb = lambda(b);
lc = lambda(c);

bx = la.*lb.*gc(:,1) - la.*lc.*gb(:,1);
by = la.*lb.*gc(:,2) - la.*lc.*gb(:,2);
bz = la.*lb.*gc(:,3) - la.*lc.*gb(:,3);

g1x = lb.*ga(:,1) + la.*gb(:,1);
g1y = lb.*ga(:,2) + la.*gb(:,2);
g1z = lb.*ga(:,3) + la.*gb(:,3);
g2x = lc.*ga(:,1) + la.*gc(:,1);
g2y = lc.*ga(:,2) + la.*gc(:,2);
g2z = lc.*ga(:,3) + la.*gc(:,3);

cx = (g1y.*gc(:,3) - g1z.*gc(:,2)) - (g2y.*gb(:,3) - g2z.*gb(:,2));
cy = (g1z.*gc(:,1) - g1x.*gc(:,3)) - (g2z.*gb(:,1) - g2x.*gb(:,3));
cz = (g1x.*gc(:,2) - g1y.*gc(:,1)) - (g2x.*gb(:,2) - g2y.*gb(:,1));
end
