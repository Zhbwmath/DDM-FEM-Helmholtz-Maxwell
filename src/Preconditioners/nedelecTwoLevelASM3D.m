function [applyPrecon, info] = nedelecTwoLevelASM3D(A_ff, edgeParts, freeEdges, P_H, opts)
% NEDELECTWOLEVELASM3D  Build two-level ASM for 3D NE_1 edge DOFs.

if nargin < 5 || isempty(opts), opts = struct(); end

A_H = P_H' * A_ff * P_H;
[coarseSolver, coarseType] = factorCoarseSPD(A_H);
[applyFine, fineInfo] = nedelecAdditiveSchwarz3D(A_ff, edgeParts, freeEdges, opts);

info = struct();
info.coarseMatrix = A_H;
info.coarseSolverType = coarseType;
info.fine = fineInfo;

    function x = applyImpl(r)
        r_H = P_H' * r;
        x_H = applyCoarse(coarseSolver, coarseType, r_H);
        x = P_H * x_H + applyFine(r);
    end

applyPrecon = @applyImpl;
end


function [solver, solverType] = factorCoarseSPD(A)
[R, p] = chol(A);
if p == 0
    solver = R;
    solverType = "chol";
else
    [L, U, P] = lu(A);
    solver = {L, U, P};
    solverType = "lu";
end
end


function x = applyCoarse(solver, solverType, r)
if solverType == "chol"
    x = solver \ (solver' \ r);
else
    x = solver{2} \ (solver{1} \ (solver{3} * r));
end
end
