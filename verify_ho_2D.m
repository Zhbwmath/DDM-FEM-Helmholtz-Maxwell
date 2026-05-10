% VERIFY_HO_2D  Convergence study for P1/P2/P3 Lagrange elements in 2D.
%
%   Solves  -\Delta u = f  on [0,1]^2  with  u = u_exact on boundary.
%   Manufactured solution:  u = sin(pi*x) * sin(pi*y)
%
%   Produces a table of L2 and H1 errors vs. mesh size h for each degree.

fprintf('========== 2D Higher-Order Convergence Study ==========\n\n');

u_exact = @(x, y) sin(pi*x) .* sin(pi*y);
f_rhs   = @(x, y) 2*pi^2 * sin(pi*x) .* sin(pi*y);

degrees = [1, 2, 3];
nRefine = 4;

fprintf('%-4s  %-8s  %-8s  %-12s  %-8s  %-12s  %-8s\n', ...
    'Deg', 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_H1', 'rateH1');
fprintf('%s\n', repmat('-', 1, 80));

for d_idx = 1:length(degrees)
    deg = degrees(d_idx);

    for k = 1:nRefine
        hk = 2^(-k-1);
        [nd, el, bd] = squaremesh([0, 1, 0, 1], hk);

        if deg > 1
            [nd_e, el_e] = extendMesh2D(nd, el, deg);
            bdNodes = getBoundaryNodesHO2D(el_e, bd, deg);
        else
            nd_e = nd;  el_e = el;
            bdNodes = getBoundaryNodesP12D(el, bd);
        end

        N_total = size(nd_e, 1);
        Ak = assembleStiffness2D(nd_e, el_e, deg);
        Mk = assembleMass2D(nd_e, el_e, deg);

        % RHS and boundary conditions
        bx = f_rhs(nd_e(:,1), nd_e(:,2));
        bk = Mk * bx;
        freeNodes = setdiff(1:N_total, bdNodes)';
        u_ex = u_exact(nd_e(:,1), nd_e(:,2));

        u_f = Ak(freeNodes, freeNodes) \ (bk(freeNodes) - Ak(freeNodes, bdNodes) * u_ex(bdNodes));
        uh = zeros(N_total, 1);
        uh(bdNodes) = u_ex(bdNodes);
        uh(freeNodes) = u_f;

        e_vec = uh - u_ex;
        errL2_k = sqrt(e_vec' * Mk * e_vec);
        errH1_k = sqrt(e_vec' * Ak * e_vec);

        if k > 1
            rL2 = log(errL2_k/errL2_prev) / log(hk/h_prev);
            rH1 = log(errH1_k/errH1_prev) / log(hk/h_prev);
            if k == 1 || d_idx == 1 && k == 1
                % header already printed
            end
            fprintf('%-4d  %-8.4f  %-8d  %-12.4e  %-8.2f  %-12.4e  %-8.2f\n', ...
                deg, hk, N_total, errL2_k, rL2, errH1_k, rH1);
        else
            fprintf('%-4d  %-8.4f  %-8d  %-12.4e  %-8s  %-12.4e  %-8s\n', ...
                deg, hk, N_total, errL2_k, '-', errH1_k, '-');
        end

        errL2_prev = errL2_k;  errH1_prev = errH1_k;  h_prev = hk;
    end
    if d_idx < length(degrees)
        fprintf('%s\n', repmat('-', 1, 80));
    end
end

fprintf('\nExpected:  P1 L2~O(h^2) H1~O(h)   P2 L2~O(h^3) H1~O(h^2)   P3 L2~O(h^4) H1~O(h^3)\n');
fprintf('========== Done ==========\n');


% ===========================================================================
function bdNode = getBoundaryNodesP12D(elem, bdFlag)
bdElem = any(bdFlag == 1, 2);
bdElemIdx = find(bdElem);
bdNode = [];
for k = 1:3
    switch k
        case 1, ev = elem(bdElemIdx, [2,3]);
        case 2, ev = elem(bdElemIdx, [3,1]);
        case 3, ev = elem(bdElemIdx, [1,2]);
    end
    isBd = bdFlag(bdElemIdx, k) == 1;
    bdNode = [bdNode; ev(isBd, :)]; %#ok<AGROW>
end
bdNode = unique(bdNode(:));
end


function bdNode = getBoundaryNodesHO2D(elem, bdFlag, degree)
% Boundary DOF indices for P2/P3 extended elements.
bdElem = any(bdFlag == 1, 2);
bdElemIdx = find(bdElem);
bdNode = [];

if degree == 2
    for k = 1:3
        isBd = bdFlag(bdElemIdx, k) == 1;
        if ~any(isBd), continue; end
        idx = bdElemIdx(isBd);
        switch k
            case 1  % edge (v2,v3): DOF columns 2,3,5
                bdNode = [bdNode; elem(idx,2); elem(idx,3); elem(idx,5)]; %#ok<AGROW>
            case 2  % edge (v3,v1): DOF columns 3,1,6
                bdNode = [bdNode; elem(idx,3); elem(idx,1); elem(idx,6)]; %#ok<AGROW>
            case 3  % edge (v1,v2): DOF columns 1,2,4
                bdNode = [bdNode; elem(idx,1); elem(idx,2); elem(idx,4)]; %#ok<AGROW>
        end
    end
else  % degree == 3
    for k = 1:3
        isBd = bdFlag(bdElemIdx, k) == 1;
        if ~any(isBd), continue; end
        idx = bdElemIdx(isBd);
        switch k
            case 1  % edge (v2,v3): cols 2,3,6,7
                bdNode = [bdNode; elem(idx,2); elem(idx,3); elem(idx,6); elem(idx,7)]; %#ok<AGROW>
            case 2  % edge (v3,v1): cols 3,1,8,9
                bdNode = [bdNode; elem(idx,3); elem(idx,1); elem(idx,8); elem(idx,9)]; %#ok<AGROW>
            case 3  % edge (v1,v2): cols 1,2,4,5
                bdNode = [bdNode; elem(idx,1); elem(idx,2); elem(idx,4); elem(idx,5)]; %#ok<AGROW>
        end
    end
end
bdNode = unique(bdNode(:));
end
