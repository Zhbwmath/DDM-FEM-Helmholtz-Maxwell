function [face, faceIdx, faceTrans] = faceMesh3D(node, elem)
% FACEMESH3D  Build global face list with NE_2 tangential trace transforms.
%
%   [face, faceIdx, faceTrans] = FACEMESH3D(node, elem)
%
%   faceTrans{t,f} is a 2-by-2 matrix T such that the two globally oriented
%   face basis traces, restricted to element t, equal the two element-local
%   face basis traces times T. The fit is done only on tangential traces,
%   since H(curl) conformity does not constrain the normal component.

NT = size(elem, 1);
faceDefs = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
nLocal = 4;

allFaces = zeros(NT*nLocal, 3);
for k = 1:nLocal
    allFaces((k-1)*NT+(1:NT), :) = elem(:, faceDefs{k});
end

[face, ~, ifa] = unique(sort(allFaces, 2), 'rows');
faceIdx = reshape(ifa, NT, nLocal);
faceTrans = cell(NT, nLocal);

for t = 1:NT
    Gphys = tetGradientsSingle(node(elem(t,:), :));

    for f = 1:nLocal
        lfIdx = faceDefs{f};
        localV = elem(t, lfIdx);
        globalV = face(faceIdx(t, f), :);

        [~, perm] = ismember(globalV, localV);
        if any(perm == 0)
            error('faceMesh3D:badFaceMap', 'Unable to map global face to local element.');
        end

        xg = node(globalV, :);
        tau1 = xg(2,:) - xg(1,:);
        tau2 = xg(3,:) - xg(1,:);

        lamVals = [1/3 1/3 1/3; 2/3 1/3 0; 0 2/3 1/3; ...
                   1/2 1/2 0; 1/4 1/2 1/4; 1/4 1/4 1/2];
        nPt = size(lamVals, 1);
        Lmat = zeros(2*nPt, 2);
        Gmat = zeros(2*nPt, 2);
        row = 0;

        for pt = 1:nPt
            lam = zeros(1, 4);
            lam(lfIdx) = lamVals(pt, :);

            L = faceBasisPair(lam, Gphys, lfIdx);
            gfIdx = lfIdx(perm);
            G = faceBasisPair(lam, Gphys, gfIdx);

            row = row + 1;
            Lmat(row, :) = [dot(L(1,:), tau1), dot(L(2,:), tau1)];
            Gmat(row, :) = [dot(G(1,:), tau1), dot(G(2,:), tau1)];

            row = row + 1;
            Lmat(row, :) = [dot(L(1,:), tau2), dot(L(2,:), tau2)];
            Gmat(row, :) = [dot(G(1,:), tau2), dot(G(2,:), tau2)];
        end

        faceTrans{t, f} = Lmat \ Gmat;
    end
end
end


function G = tetGradientsSingle(v)
e12 = v(2,:) - v(1,:);
e13 = v(3,:) - v(1,:);
e14 = v(4,:) - v(1,:);
J = [e12; e13; e14]';
invJ = inv(J);
g2 = invJ(1,:);
g3 = invJ(2,:);
g4 = invJ(3,:);
g1 = -(g2 + g3 + g4);
G = {g1, g2, g3, g4};
end


function B = faceBasisPair(lambda, Gphys, idx)
i = idx(1);
j = idx(2);
k = idx(3);

B = zeros(2, 3);
B(1, :) = lambda(i) * lambda(j) * Gphys{k} - lambda(i) * lambda(k) * Gphys{j};
B(2, :) = lambda(j) * lambda(k) * Gphys{i} - lambda(j) * lambda(i) * Gphys{k};
end
