function [applyPrecon, info] = nedelecAdditiveSchwarz3D(A_ff, edgeParts, freeEdges, opts)
% NEDELECADDITIVESCHWARZ3D  Build one-level ASM for 3D NE_1 edge DOFs.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'solver'), opts.solver = 'chol'; end

Nf = numel(freeEdges);
global2reduced = zeros(max(freeEdges), 1);
global2reduced(freeEdges) = (1:Nf).';

nSub = numel(edgeParts);
locSolvers = cell(nSub, 1);
locRedIdx = cell(nSub, 1);
solverType = strings(nSub, 1);

for s = 1:nSub
    interior = edgeParts(s).interiorEdgeIdx(:);
    freeInterior = intersect(interior, freeEdges, 'stable');
    if isempty(freeInterior)
        locRedIdx{s} = [];
        locSolvers{s} = [];
        solverType(s) = "empty";
        continue;
    end

    redPos = global2reduced(freeInterior);
    A_loc = A_ff(redPos, redPos);
    [locSolvers{s}, solverType(s)] = factorLocalSPD(A_loc, opts.solver);
    locRedIdx{s} = redPos;
end

info = struct();
info.localReducedDofs = locRedIdx;
info.localSolverType = solverType;
info.freeEdges = freeEdges;

    function x = applyImpl(r)
        x = zeros(Nf, size(r, 2));
        for ss = 1:nSub
            redPos = locRedIdx{ss};
            if isempty(redPos), continue; end
            x(redPos, :) = x(redPos, :) + applyLocal(locSolvers{ss}, solverType(ss), r(redPos, :));
        end
    end

applyPrecon = @applyImpl;
end


function [solver, solverType] = factorLocalSPD(A, preferred)
if strcmpi(preferred, 'chol')
    [R, p] = chol(A);
    if p == 0
        solver = R;
        solverType = "chol";
        return;
    end
end

[L, U, P] = lu(A);
solver = {L, U, P};
solverType = "lu";
end


function x = applyLocal(solver, solverType, r)
if solverType == "chol"
    x = solver \ (solver' \ r);
else
    x = solver{2} \ (solver{1} \ (solver{3} * r));
end
end
