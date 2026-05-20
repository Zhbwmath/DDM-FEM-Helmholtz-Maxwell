function applyPrecon = additiveSchwarz(A_ff, partitions, freeNodes)
% ADDITIVESCHWARZ  Build ASM preconditioner with Dirichlet inner BC.
%
%   M^{-1} = sum_i R_i^T A_i^{-1} R_i
%
%   Each subdomain Ω_i uses homogeneous Dirichlet on ∂Ω_i \ ∂Ω.
%   R_i restricts global free DOFs to interior free DOFs of Ω_i.
%   V_{h,i} = {v ∈ V_h|_{Ω_i} : v = 0 on ∂Ω_i}
%
%   applyPrecon = ADDITIVESCHWARZ(A_ff, partitions, freeNodes)
%
%   Input:
%     A_ff       - reduced global stiffness A(freeNodes, freeNodes)
%     partitions - from partitionMesh with 'overlap', delta (has .interiorNodeIdx)
%     freeNodes  - list of global free DOF indices
%   Output:
%     applyPrecon - function handle: x = applyPrecon(r), r,x of length Nf

Nf = length(freeNodes);
global2reduced = zeros(max(freeNodes), 1);
global2reduced(freeNodes) = (1:Nf)';

nSub = length(partitions);
locSolvers = cell(nSub, 1);
locRedIdx = cell(nSub, 1);  % positions in the reduced global vector

for s = 1:nSub
    % Interior nodes of Ω_i that are also global free DOFs
    interior = partitions(s).interiorNodeIdx;
    freeInterior = intersect(interior, freeNodes);

    if isempty(freeInterior)
        locSolvers{s} = [];
        locRedIdx{s} = [];
        continue;
    end

    % Their positions in the reduced global vector
    redPos = global2reduced(freeInterior);

    % Extract subdomain matrix and factor
    A_loc = A_ff(redPos, redPos);
    try
        locSolvers{s} = chol(A_loc);  % Cholesky for SPD
    catch
        [L, U, P] = lu(A_loc);
        locSolvers{s} = {L, U, P};
    end

    locRedIdx{s} = redPos;
end

    function x = applyImpl(r)
        x = zeros(Nf, 1);
        for s = 1:nSub
            redPos = locRedIdx{s};
            if isempty(redPos), continue; end

            r_loc = r(redPos);
            solver = locSolvers{s};

            if isempty(solver)
                continue;
            elseif ismatrix(solver) && size(solver, 1) == size(solver, 2)
                x_loc = solver \ (solver' \ r_loc);
            else
                x_loc = solver{2} \ (solver{1} \ (solver{3} * r_loc));
            end

            x(redPos) = x(redPos) + x_loc;
        end
    end

applyPrecon = @applyImpl;
end
