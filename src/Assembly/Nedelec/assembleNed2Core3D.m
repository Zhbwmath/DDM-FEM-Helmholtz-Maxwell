function A = assembleNed2Core3D(node, elem, op)
% ASSEMBLENED2CORE3D  Assemble conforming 3D NE_2 mass or curl-curl matrix.

[gIdx, trans, Ntot] = ned2Dof3D(node, elem);
NT = size(elem, 1);
nLocal = 20;

[lambda, weight] = quadtet(4);
nQuad = length(weight);
nEntries = NT * nLocal * nLocal * nQuad;

ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    [Bx, By, Bz, Cx, Cy, Cz, volume] = ned2TransformedBasis3D(node, elem, lambda(q,:), trans);
    w = 6 * weight(q) * volume;

    for p = 1:nLocal
        gp = gIdx(:, p);
        for r = 1:nLocal
            gr = gIdx(:, r);
            switch lower(op)
                case 'mass'
                    val = Bx(:,p).*Bx(:,r) + By(:,p).*By(:,r) + Bz(:,p).*Bz(:,r);
                case 'curl'
                    val = Cx(:,p).*Cx(:,r) + Cy(:,p).*Cy(:,r) + Cz(:,p).*Cz(:,r);
                otherwise
                    error('assembleNed2Core3D:badOp', 'Unknown operator: %s', op);
            end
            nxt = idx + 1;
            idx = idx + NT;
            ii(nxt:idx) = gp;
            jj(nxt:idx) = gr;
            ss(nxt:idx) = w .* val;
        end
    end
end

A = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Ntot, Ntot);
end
