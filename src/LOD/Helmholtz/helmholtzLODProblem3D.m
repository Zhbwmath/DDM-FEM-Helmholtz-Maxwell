function problem = helmholtzLODProblem3D(nodeH, elemH, ~, nodeh, elemh, bdh, k, f, g)
% HELMHOLTZLODPROBLEM3D  Matrix-callback data for 3D Helmholtz LOD.

if nargin < 8, f = []; end
if nargin < 9, g = []; end

problem = struct();
problem.bdFlagFine = bdh;
problem.transfer = @() prolongateNestedP1(nodeH, elemH, nodeh);
problem.interpolation = @() weightedClementP1(nodeh, elemh, nodeH, elemH);
problem.constraints.patch = @(Q, patch, T, opts) lodClementConstraints(Q, patch, T, opts);
problem.form.global = @globalForm;
problem.form.patch = @patchForm;
problem.form.elementRhs = @elementRhs;
problem.form.elementRhsAdjoint = @elementRhsAdjoint;

    function [A, b] = globalForm()
        [A, b] = assembleHelmholtz3D(nodeh, elemh, bdh, k, f, g, 1);
    end

    function A = patchForm(patch, T)
        A = assembleHelmholtz3D(patch.localNode{T}, patch.localElem{T}, ...
            patch.localBdFlag{T}, k, [], [], 1);
    end

    function R = elementRhs(Tcoarse, targetDof, patch, Tpatch, P)
        R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, P, false);
    end

    function R = elementRhsAdjoint(Tcoarse, targetDof, patch, Tpatch, P)
        R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, P, true);
    end

    function R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, P, adjoint)
        fineIds = patch.targetFineElemIds{Tcoarse};
        local2global = patch.local2global{Tpatch};
        [isInPatch, localElem] = ismember(elemh(fineIds, :), local2global);
        if any(~isInPatch(:))
            error('helmholtzLODProblem3D:patchMap', ...
                'Target fine element is not contained in its patch.');
        end

        Aelem = assembleHelmholtz3D(patch.localNode{Tpatch}, localElem, ...
            bdh(fineIds, :), k, [], [], 1);
        V = P(local2global, targetDof);
        if adjoint
            R = Aelem' * V;
        else
            R = Aelem * V;
        end
    end
end
