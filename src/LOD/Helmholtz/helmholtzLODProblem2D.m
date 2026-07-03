function problem = helmholtzLODProblem2D(nodeH, elemH, ~, nodeh, elemh, bdh, k, f, g)
% HELMHOLTZLODPROBLEM2D  Matrix-callback data for Helmholtz LOD.

if nargin < 8, f = []; end
if nargin < 9, g = []; end
kInput = helmholtzLODInput(k);

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
        [A, b] = assembleHelmholtz2D(nodeh, elemh, bdh, kInput, f, g, 1);
    end

    function A = patchForm(patch, T)
        sub = lodGetPatchSubmesh(patch, T);
        A = assembleHelmholtz2D(sub.localNode, sub.localElem, ...
            sub.localBdFlag, kInput, [], [], 1);
    end

    function R = elementRhs(Tcoarse, targetDof, patch, Tpatch, P)
        R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, P, false);
    end

    function R = elementRhsAdjoint(Tcoarse, targetDof, patch, Tpatch, P)
        R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, P, true);
    end

    function R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, P, adjoint)
        fineIds = patch.targetFineElemIds{Tcoarse};
        sub = lodGetPatchSubmesh(patch, Tpatch);
        local2global = sub.local2global;
        [isInPatch, localElem] = ismember(elemh(fineIds, :), local2global);
        if any(~isInPatch(:))
            error('helmholtzLODProblem2D:patchMap', ...
                'Target fine element is not contained in its patch.');
        end

        Aelem = assembleHelmholtz2D(sub.localNode, localElem, ...
            bdh(fineIds, :), kInput, [], [], 1);
        V = P(local2global, targetDof);
        if adjoint
            R = Aelem' * V;
        else
            R = Aelem * V;
        end
    end
end


function input = helmholtzLODInput(k)
if isnumeric(k) && isscalar(k)
    input = k;
else
    input = normalizeHelmholtzPDE(k);
end
end
