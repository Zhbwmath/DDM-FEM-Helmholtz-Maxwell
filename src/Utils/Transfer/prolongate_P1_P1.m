function P = prolongate_P1_P1(coarseNode, coarseElem, fineNode)
% PROLONGATE_P1_P1  P1 prolongation from coarse mesh to fine mesh.
%
%   P = PROLONGATE_P1_P1(coarseNode, coarseElem, fineNode)
%
%   For each fine node x_k, finds the coarse element containing it and
%   evaluates the P1 (linear) basis functions at x_k.
%   P is sparse, size N_fine × N_coarse.

Nc = size(coarseNode, 1);
Nf = size(fineNode, 1);
NTc = size(coarseElem, 1);

% Preallocate triplets (at most 3*Nf entries for P1)
ii = zeros(3 * Nf, 1);
jj = zeros(3 * Nf, 1);
ss = zeros(3 * Nf, 1);
idx = 0;

for k = 1:Nf
    xk = fineNode(k, :);

    % Find which coarse element contains xk (brute force for simplicity)
    found = false;
    for e = 1:NTc
        v = coarseElem(e, :);
        x1 = coarseNode(v(1), :);  x2 = coarseNode(v(2), :);  x3 = coarseNode(v(3), :);

        % Barycentric coordinates
        denom = (x2(2)-x3(2))*(x1(1)-x3(1)) + (x3(1)-x2(1))*(x1(2)-x3(2));
        lam1 = ((x2(2)-x3(2))*(xk(1)-x3(1)) + (x3(1)-x2(1))*(xk(2)-x3(2))) / denom;
        lam2 = ((x3(2)-x1(2))*(xk(1)-x3(1)) + (x1(1)-x3(1))*(xk(2)-x3(2))) / denom;
        lam3 = 1 - lam1 - lam2;

        if lam1 >= -1e-12 && lam2 >= -1e-12 && lam3 >= -1e-12
            % xk is in this element
            idx = idx + 1; ii(idx) = k; jj(idx) = v(1); ss(idx) = lam1;
            idx = idx + 1; ii(idx) = k; jj(idx) = v(2); ss(idx) = lam2;
            idx = idx + 1; ii(idx) = k; jj(idx) = v(3); ss(idx) = lam3;
            found = true;
            break;
        end
    end
    if ~found
        warning('Fine node %d at (%.4f,%.4f) not found in any coarse element', k, xk);
    end
end

P = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Nf, Nc);
end
