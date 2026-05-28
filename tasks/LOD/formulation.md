Created: 2026-05-26
Updated: 2026-05-27
Verification entry point: `verify/verify_lod_corrector_orthogonality.m`
Main utilities: `buildLOD`, `lodBuildPatches`, `lodClementConstraints`, `lodSolveConstrainedSaddle`

# Abstract LOD Framework

This note records the reusable localized orthogonal decomposition framework used by the Helmholtz reproduction in `tasks/LOD/LOD_Helmholtz`.

Let \(V_h\) be a fine conforming finite element space and \(V_H\subset V_h\) a coarse P1 space. A quasi-interpolation \(I_H:V_h\to V_H\) defines the constrained fine-scale space
\[
W_h=\ker I_H.
\]
For a problem-adapted sesquilinear form \(a(\cdot,\cdot)\), the localized primal corrector for a coarse element \(T\) and patch \(P=\Omega_{T,\ell}\) is computed from
\[
a_P(q_{T,y},w)=a_T(\phi_y,w)\quad\forall w\in W_h(P).
\]
The adjoint corrector uses
\[
a_P(w,q^*_{T,y})=a_T(w,\phi_y)\quad\forall w\in W_h(P).
\]

The implementation does not construct a basis of \(W_h(P)\) for production runs. Instead, homogeneous constraints are represented by a local matrix \(C_P\), with
\[
C_P^* q=0.
\]
The local corrector block is solved as the saddle system
```text
[ A_P   C_P ] [ Q_T      ] = [ R_T ]
[ C_P'   0  ] [ Lambda_T ]   [  0  ]
```
where `v' * A_P * u = a_P(u,v)`. Explicit nullspace bases are allowed only in small verification scripts.

Corrected trial and test bases are
\[
\psi_z=\phi_z-\sum_{T\ni z}q_{T,z},\qquad
\psi^*_z=\phi_z-\sum_{T\ni z}q^*_{T,z}.
\]
The Petrov-Galerkin coarse matrix is assembled by
\[
(A_H)_{ij}=a(\psi_j,\psi^*_i).
\]

The patch layer now covers nested P1 triangle and tetrahedron meshes. In both dimensions, artificial patch-boundary nodes are fixed out of the local corrector problem, while physical boundary nodes remain free so the local form can retain the impedance boundary contribution.
