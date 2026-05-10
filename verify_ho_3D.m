% VERIFY_HO_3D  Convergence study for P1/P2/P3 in 3D.
%   Solves -Delta u = f on [0,1]^3 with u = u_exact on boundary.
%   u = sin(pi*x)*sin(pi*y)*sin(pi*z),  f = 3*pi^2*u.

fprintf('========== 3D Higher-Order Convergence Study ==========\n\n');

u_exact = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z);
f_rhs   = @(x,y,z) 3*pi^2*u_exact(x,y,z);

degrees = [1, 2, 3];
nRefine = 3;

fprintf('%-4s  %-8s  %-8s  %-12s  %-8s  %-12s  %-8s\n', ...
    'Deg', 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_H1', 'rateH1');
fprintf('%s\n', repmat('-', 1, 80));

for d_idx = 1:length(degrees)
    deg = degrees(d_idx);
    for k = 1:nRefine
        hk = 2^(-k-1);
        [nd, el, bd] = cubemesh([0,1,0,1,0,1], hk);
        if deg > 1
            [nd_e, el_e] = extendMesh3D(nd, el, deg);
            bdNodes = getBd3D(el_e, bd, deg);
        else
            nd_e = nd;  el_e = el;
            bdNodes = getBd3D(el, bd, 1);
        end
        N = size(nd_e, 1);
        Ak = assembleStiffness3D(nd_e, el_e, deg);
        Mk = assembleMass3D(nd_e, el_e, deg);
        bk = Mk * f_rhs(nd_e(:,1), nd_e(:,2), nd_e(:,3));
        free = setdiff(1:N, bdNodes)';
        u_ex = u_exact(nd_e(:,1), nd_e(:,2), nd_e(:,3));
        u_f = Ak(free,free) \ (bk(free) - Ak(free,bdNodes)*u_ex(bdNodes));
        uh = zeros(N,1); uh(bdNodes)=u_ex(bdNodes); uh(free)=u_f;
        e = uh - u_ex;
        eL2 = sqrt(e'*Mk*e);  eH1 = sqrt(e'*Ak*e);
        if k>1
            rL2 = log(eL2/eL2p)/log(hk/hp);
            rH1 = log(eH1/eH1p)/log(hk/hp);
            fprintf('%-4d  %-8.4f  %-8d  %-12.4e  %-8.2f  %-12.4e  %-8.2f\n', ...
                deg, hk, N, eL2, rL2, eH1, rH1);
        else
            fprintf('%-4d  %-8.4f  %-8d  %-12.4e  %-8s  %-12.4e  %-8s\n', ...
                deg, hk, N, eL2, '-', eH1, '-');
        end
        eL2p=eL2; eH1p=eH1; hp=hk;
    end
    if d_idx<length(degrees), fprintf('%s\n', repmat('-',1,80)); end
end
fprintf('\nExpected: P1 L2~O(h^2) H1~O(h)  P2 L2~O(h^3) H1~O(h^2)  P3 L2~O(h^4) H1~O(h^3)\n');
fprintf('========== Done ==========\n');

function bdNode = getBd3D(elem, bdFlag, degree)
fv = {[2,3,4],[1,4,3],[1,2,4],[1,3,2]};
bdNode = [];
if degree == 1
    for f=1:4
        isF = bdFlag(:,f)==1;
        if ~any(isF), continue; end
        bdNode = [bdNode; elem(isF,fv{f}(1)); elem(isF,fv{f}(2)); elem(isF,fv{f}(3))]; %#ok<AGROW>
    end
elseif degree == 2
    fd = {[2,3,4,8,10,9], [1,3,4,6,10,7], [1,2,4,5,9,7], [1,2,3,5,8,6]};
    for f=1:4
        isF = bdFlag(:,f)==1;
        if ~any(isF), continue; end
        bdNode = [bdNode; reshape(elem(isF,fd{f}),[],1)]; %#ok<AGROW>
    end
else
    fd = {[2,3,4,11,12,15,16,14,13,20], [1,3,4,7,8,15,16,10,9,19], ...
          [1,2,4,5,6,13,14,9,10,18],   [1,2,3,5,6,11,12,8,7,17]};
    for f=1:4
        isF = bdFlag(:,f)==1;
        if ~any(isF), continue; end
        bdNode = [bdNode; reshape(elem(isF,fd{f}),[],1)]; %#ok<AGROW>
    end
end
bdNode = unique(bdNode(:));
end
