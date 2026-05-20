function applyPrecon = twoLevelASM(A_ff, parts, freeNodes, P_H)
% TWOLEVELASM  Two-level additive Schwarz preconditioner.
%
%   M^{-1} = P_H A_H^{-1} R_H + Σ_i R_i^T A_i^{-1} R_i
%
%   applyPrecon = TWOLEVELASM(A_ff, parts, freeNodes, P_H)
%
%   P_H: prolongation from coarse to fine (Nf × Nc).
%        Assumed to be constructed on freeNodes only.

Nf = size(A_ff, 1);

% Coarse matrix (Galerkin)
R_H = P_H';
A_H = R_H * A_ff * P_H;

% Factor coarse matrix
try
    L_H = chol(A_H);
    coarseSolver = 'chol';
catch
    [L_H, U_H, P_H_lu] = lu(A_H);
    coarseSolver = 'lu';
end

% Fine-level ASM
applyFine = additiveSchwarz(A_ff, parts, freeNodes);

    function x = applyImpl(r)
        % Coarse correction
        r_H = R_H * r;
        if strcmp(coarseSolver, 'chol')
            x_H = L_H \ (L_H' \ r_H);
        else
            x_H = U_H \ (L_H \ (P_H_lu * r_H));
        end
        x = P_H * x_H;

        % Add fine-level correction
        x = x + applyFine(r);
    end

applyPrecon = @applyImpl;
end
