% VERIFY_NED2_TRACE3D  Check NE_2 face tangential trace conformity in 3D.

fprintf('========== 3D NE_2 Face Trace Conformity ==========\n\n');

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

maxJump = checkNed2Trace3D(0.5);
assert(maxJump < 1e-10, ...
    'NE_2 3D face tangential trace jump is too large: %.3e', maxJump);

fprintf('Tangential face traces are continuous to roundoff.\n');
fprintf('========== Done ==========\n');


function maxJump = checkNed2Trace3D(h)
[node, elem] = cubemesh([0,1,0,1,0,1], h);
[face, faceIdx] = faceMesh3D(node, elem);
[~, trans] = ned2Dof3D(node, elem);

NF = size(face, 1);
owners = cell(NF, 1);
for t = 1:size(elem, 1)
    for f = 1:4
        fid = faceIdx(t, f);
        owners{fid} = [owners{fid}; t, f]; %#ok<AGROW>
    end
end

faceDefs = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
lambdaFace = [1/3, 1/3, 1/3];
maxJump = 0;

for fid = 1:NF
    if size(owners{fid}, 1) ~= 2
        continue
    end

    xg = node(face(fid, :), :);
    tau1 = xg(2,:) - xg(1,:);
    tau2 = xg(3,:) - xg(1,:);

    for jdof = 1:2
        tr = zeros(2, 2);
        for side = 1:2
            t = owners{fid}(side, 1);
            f = owners{fid}(side, 2);
            lf = faceDefs{f};
            localFaceVertices = elem(t, lf);
            [~, perm] = ismember(face(fid, :), localFaceVertices);

            lambda = zeros(1, 4);
            lambda(lf(perm)) = lambdaFace;

            [Bx, By, Bz] = ned2TransformedBasis3D(node, elem(t, :), lambda, trans(t, :, :));
            col = 12 + 2*(f-1) + jdof;
            v = [Bx(1, col), By(1, col), Bz(1, col)];
            tr(side, :) = [dot(v, tau1), dot(v, tau2)];
        end
        maxJump = max(maxJump, norm(tr(1,:) - tr(2,:), inf));
    end
end

fprintf('NE2 3D max face tangential trace jump: %.3e\n', maxJump);
end
