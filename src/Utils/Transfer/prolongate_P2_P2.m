function P = prolongate_P2_P2(coarseNode, coarseElem, fineNode)
% PROLONGATE_P2_P2  P2 prolongation from coarse mesh to fine mesh.
%
%   P = PROLONGATE_P2_P2(coarseNode, coarseElem, fineNode)
%
%   coarseNode, fineNode: P2 node coordinates (already extended).
%   coarseElem: P2 element connectivity (NTc × 6).

Nc = size(coarseNode, 1);
Nf = size(fineNode, 1);
NTc = size(coarseElem, 1);

ii = zeros(6 * Nf, 1);  % P2 has up to 6 DOFs per element
jj = zeros(6 * Nf, 1);
ss = zeros(6 * Nf, 1);
idx = 0;

% Precompute P2 basis at fine node positions (reuse lagrange2D)
for k = 1:Nf
    xk = fineNode(k, :);

    found = false;
    for e = 1:NTc
        v = coarseElem(e, 1:3);  % vertex indices only for barycentric check
        x1 = coarseNode(v(1), :);  x2 = coarseNode(v(2), :);  x3 = coarseNode(v(3), :);

        % Barycentric coordinates of xk in this coarse triangle
        denom = (x2(2)-x3(2))*(x1(1)-x3(1)) + (x3(1)-x2(1))*(x1(2)-x3(2));
        lam1 = ((x2(2)-x3(2))*(xk(1)-x3(1)) + (x3(1)-x2(1))*(xk(2)-x3(2))) / denom;
        lam2 = ((x3(2)-x1(2))*(xk(1)-x3(1)) + (x1(1)-x3(1))*(xk(2)-x3(2))) / denom;
        lam3 = 1 - lam1 - lam2;

        if lam1 >= -1e-12 && lam2 >= -1e-12 && lam3 >= -1e-12
            % Evaluate P2 basis at barycentric coordinates
            lambda = [lam1, lam2, lam3];
            phi = lagrange2D(2, lambda);  % 1 × 6 P2 basis values

            for j = 1:6
                idx = idx + 1;
                ii(idx) = k;
                jj(idx) = coarseElem(e, j);
                ss(idx) = phi(j);
            end
            found = true;
            break;
        end
    end
    if ~found
        warning('Fine P2 node %d at (%.4f,%.4f) not found in any coarse P2 element', k, xk);
    end
end

P = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Nf, Nc);
end
