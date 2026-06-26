Reproduction target: Hu--Li, *A novel coarse space applying to the weighted Schwarz method for Helmholtz equations*, Tables 1--3 and algorithm design.
Created: 2026-06-10
Updated: 2026-06-26
Verification entry point: `verify/verify_hl25_coarse_spaces.m`; `verify/verify_hl25_tables123.m`; `verify/verify_hl25_lxzz_cross_study.m`; `verify/verify_hl25_lxzz_hybrid_medium.m`; `verify/verify_hl25_full_sweep.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `partitionMesh2D`, `linearPartitionOfUnity2D`, `assembleHelmholtz2D`, MATLAB `gmres`

# Hu--Li Helmholtz-Harmonic Coarse Spaces

## Model Problem

On $\Omega=(0,1)^2$, solve

$$
-\Delta u-(\kappa^2+i\varepsilon)u=f
\quad\text{in }\Omega,
$$

with impedance boundary condition

$$
\partial_nu-i\kappa u=g
\quad\text{on }\partial\Omega.
$$

The experiments use the plane wave

$$
u(x)=\exp(i\kappa d\cdot x),
\qquad
d=(1/\sqrt2,1/\sqrt2),
$$

with $f=0$ and its impedance boundary data. The paper compares
$\varepsilon=0$ and $\varepsilon=\kappa$.

For a P1 or P2 nodal finite element space $V_h$, the global matrix is

$$
A_\varepsilon
=K-M_{\kappa^2+i\varepsilon}-iM_{\kappa,\partial\Omega}.
$$

The positive energy matrix is

$$
D_\kappa=K+M_{|\kappa|^2}.
$$

For variable $\kappa(x)$, all mass and boundary terms are assembled with
quadrature. A scalar `kappaRef` is additionally required for paper-based mesh,
threshold, and economic-space scaling.

## Subdomains And Weights

The domain is decomposed into overlapping rectangular subdomains
$\{\Omega_\ell\}_{\ell=1}^{N_s}$. The base subdomain width is

$$
d\approx\kappa^{-\beta}.
$$

The reproduction uses the generous-overlap setting

$$
\delta\approx d/4,
$$

while the implementation also permits minimal overlap $\delta=h$. The fine
mesh, subdomain boundaries, and overlap boundaries are aligned:

$$
d=k_dh,\qquad \delta=k_\delta h
$$

for integers $k_d$ and $k_\delta$.

Let $\chi_\ell$ be the normalized nodal partition-of-unity weight and let
$X_\ell$ be its diagonal matrix on the local degrees of freedom. The global
weighted extension is

$$
\widetilde E_\ell=E_\ell X_\ell.
$$

## Local Helmholtz-Harmonic Space

On each subdomain, assemble the full-impedance local matrix

$$
A_{\varepsilon,\ell}
=K_\ell-M_{\kappa^2+i\varepsilon,\ell}
-iM_{\kappa,\partial\Omega_\ell}
$$

and the local energy matrix

$$
D_{\kappa,\ell}
=K_\ell+M_{|\kappa|^2,\ell}.
$$

Let $E_{\partial,\ell}$ inject local boundary degrees of freedom and let
$M_{\partial,\ell}$ be the unweighted boundary mass matrix. The impedance-data
Helmholtz-harmonic basis matrix is

$$
H_\ell
=A_{\varepsilon,\ell}^{-1}
M_{\partial,\ell}E_{\partial,\ell}.
$$

Thus every column satisfies

$$
A_{\varepsilon,\ell}h
=M_{\partial,\ell}E_{\partial,\ell}\lambda.
$$

The implementation factors $A_{\varepsilon,\ell}$ once and reuses that factor
for harmonic extensions and the native weighted Schwarz action.

## Spectral Coarse Space

The Hu--Li local generalized eigenproblem becomes

$$
H_\ell^*X_\ell D_{\kappa,\ell}X_\ell H_\ell c
=\lambda H_\ell^*D_{\kappa,\ell}H_\ell c.
$$

Retain every eigenvector satisfying

$$
\lambda\geq\rho^2.
$$

The corresponding global coarse columns are

$$
E_\ell X_\ell H_\ell c.
$$

Small trace problems use dense Hermitian eigensolves. Larger trace problems use
an adaptive largest-eigenvalue `eigs` solve, increasing the requested dimension
until the first computed eigenvalue is below $\rho^2$.

The theoretical Table 1 choices are

$$
\rho=\kappa^{\beta-2-1/(4p)}
\quad\text{for }\varepsilon=\kappa,
$$

and

$$
\rho=\kappa^{-1}
\quad\text{for }\varepsilon=0.
$$

Table 2 uses

$$
\rho=\frac12\kappa^{(\beta-1)/2}.
$$

## Economic Coarse Space

Parameterize each axis-aligned rectangular $\partial\Omega_\ell$ periodically
by normalized arclength $t\in[0,1)$. Let $V_\nu([0,1])$ be the periodic
quadratic finite element space on $\nu$ uniform boundary elements. It has
$2\nu$ periodic nodal basis functions.

Evaluate those basis functions at the fine boundary degrees of freedom to
obtain a trace matrix $T_{\nu,\ell}$. The local economic basis is

$$
H_\ell T_{\nu,\ell},
$$

and its weighted global contribution is

$$
E_\ell X_\ell H_\ell T_{\nu,\ell}.
$$

The implementation does not form every column of $H_\ell$ for the economic
space. It applies the same operator more efficiently by solving

$$
A_{\varepsilon,\ell}Y_{\nu,\ell}
=M_{\partial,\ell}E_{\partial,\ell}T_{\nu,\ell}
$$

directly for the $2\nu$ required columns. This is algebraically identical to
$Y_{\nu,\ell}=H_\ell T_{\nu,\ell}$ and avoids extending unused boundary data.

Table 3 uses

$$
h_\nu=\kappa^{\beta-1},
\qquad
\nu=\operatorname{round}(\kappa^{1-\beta}).
$$

## Global Coarse Solver

Concatenate all weighted local basis columns into $Z_{\rm raw}$. The columns
are normalized, and a rank-revealing QR is applied to the small Gram matrix

$$
G=\widehat Z_{\rm raw}^*\widehat Z_{\rm raw}.
$$

The selected original columns remain sparse. This removes numerical
dependencies without the expensive tall QR of the full fine-by-coarse matrix
and without changing the represented coarse space:

$$
\operatorname{range}(Z)=\operatorname{range}(Z_{\rm raw}).
$$

The coarse matrix and residual inverse are

$$
A_0=Z^*A_\varepsilon Z,
\qquad
M_0^{-1}=ZA_0^{-1}Z^*.
$$

The object returned in `method.coarseSpace` has `trial=test=Z`, `AH=A0`, and
coarse solve/adjoint-solve handles. It can be passed directly through
`opts.coarseSpace` to `twoLevelHybridSchwarzHelmholtz2D`.

## Native Hu--Li Hybrid

The one-level weighted additive Schwarz inverse is

$$
M_{\rm WASI}^{-1}
=\sum_{\ell=1}^{N_s}
E_\ell X_\ell A_{\varepsilon,\ell}^{-1}E_\ell^*.
$$

The residual-level native hybrid preconditioner is

$$
B_{\rm HL}^{-1}
=(I-M_0^{-1}A_\varepsilon)M_{\rm WASI}^{-1}+M_0^{-1}.
$$

For a residual $r$, the implementation computes

$$
z_0=M_0^{-1}r,\qquad
z_w=M_{\rm WASI}^{-1}r,\qquad
z=z_0+z_w-M_0^{-1}A_\varepsilon z_w.
$$

## LXZZ25 Integration

The same `fineSpace` and `coarseSpace` objects are injected into the existing
LXZZ25 wrapper. No change is made to its hybrid algebra.

For the energy-adjoint variant, the factorization of $D_\kappa$ is cached in
the injected fine-space object and shared by both LXZZ local variants and both
Hu--Li coarse constructions on the same mesh. An optional precomputed matrix
$D_\kappa^{-1}A_\varepsilon^*Z$ is supported and is algebraically exact, but
it is disabled by default because profiling showed that its many-right-hand-
side setup cost exceeds the saved Krylov-application time for the tested
iteration counts.

The cross-study compares:

1. Native Hu--Li weighted hybrid.
2. LXZZ25 Dirichlet local hybrid.
3. LXZZ25 impedance local hybrid.

### Resolved LXZZ Operator Form

This subsection uses LXZZ notation: $P_0$ is the coarse-to-fine injection, not
a hybrid complement. The injected Hu--Li coarse space supplies a trial
injection $P_0$ and a test injection $P_0^\star$; in the current Hu--Li
injection $P_0^\star=P_0=Z$. Let $A=A_\varepsilon$ be the fine Helmholtz
matrix and let $D=D_\kappa$ be the Hermitian positive energy matrix. Define

$$
A_0=(P_0^\star)^*AP_0,\qquad
M_0^{-1}=P_0A_0^{-1}(P_0^\star)^*,
$$

and the function-level coarse Ritz operator

$$
Q_0=M_0^{-1}A
=P_0A_0^{-1}(P_0^\star)^*A.
$$

The default LXZZ cross-study uses the $D_\kappa$-adjoint of $Q_0$. With
$\langle u,v\rangle_D=v^*Du$,

$$
Q_0^{T_D}=D^{-1}Q_0^*D
=D^{-1}A^*P_0^\star A_0^{-*}P_0^*D.
$$

For the current Hu--Li injection $P_0^\star=P_0=Z$, this is

$$
Q_0^{T_D}=D^{-1}A^*Z A_0^{-*}Z^*D.
$$

The LXZZ operators are twice-hybrid function-level operators. If

$$
Q_D=M_D^{-1}A,\qquad Q_I=M_I^{-1}A,
$$

then the Dirichlet-local and impedance-local variants are

$$
Q_m^{(1)}
=Q_0+(I-Q_0)^{T_D}Q_D(I-Q_0),
$$

and

$$
Q_m^{(2)}
=Q_0+(I-Q_0)^{T_D}Q_I(I-Q_0).
$$

This is distinct from the Hu--Li native one-hybrid residual formula. In code,
`apply(v)` exposes $Q_m^{(i)}v$.

For left-preconditioned GMRES, the implementation applies the residual
preconditioner $B_i^{-1}$, not $Q_m^{(i)}$ directly. The equivalent residual
forms are

$$
B_1^{-1}
=M_0^{-1}+(I-Q_0)^{T_D}M_D^{-1}(I-AM_0^{-1}),
$$

and

$$
B_2^{-1}
=M_0^{-1}+(I-Q_0)^{T_D}M_I^{-1}(I-AM_0^{-1}).
$$

They satisfy

$$
B_i^{-1}A=Q_m^{(i)}
$$

up to roundoff; this is the identity checked in code as
`apply(v) == applyResidual(A*v)`.

For the Dirichlet LXZZ variant, let $\Omega_j^{\rm LXZZ}$ be one local
subdomain of the LXZZ partition, $I_j$ the local degrees of freedom touched by
that patch, and $F_j\subset I_j$ the active interior degrees of freedom after
imposing homogeneous Dirichlet conditions on the artificial boundary. Let
$R_j:\mathbb{C}^N\to\mathbb{C}^{|F_j|}$ restrict a global vector to $F_j$.
Let $A_j^D$ be the local Helmholtz matrix assembled on
$\Omega_j^{\rm LXZZ}$ and restricted to $F_j$. Then

$$
M_{{\rm loc},D}^{-1}
=\sum_j R_j^*(A_j^D)^{-1}R_j.
$$

For the impedance LXZZ variant, all patch degrees of freedom $I_j$ are used.
Let $C_j$ be the coercive local impedance matrix

$$
C_j=K_j+\kappa^2M_j-\mathrm{i}\kappa M_{\partial\Omega_j},
$$

with the corresponding variable-coefficient matrices used when $k(x)$ is
provided. Let $W_j:\mathbb{C}^{|I_j|}\to\mathbb{C}^N$ be the weighted
extension whose nodal weights are the local hat weights normalized by the sum
of all patch hat weights at the same node. Then

$$
M_{{\rm loc},I}^{-1}
=\sum_j W_j C_j^{-1}W_j^*.
$$

Thus the two concrete preconditioners used in the cross-study are obtained by
substituting $M_{{\rm loc},D}^{-1}$ or $M_{{\rm loc},I}^{-1}$ for
$M_D^{-1}$ or $M_I^{-1}$ in the twice-hybrid formulas above.

Each comparison uses the same $A_\varepsilon$, fine mesh, and Hu--Li coarse
basis. In runs where the goal is to study the Hu--Li paper decomposition, the
Hu--Li coarse basis keeps $d\approx\kappa^{-\beta}$ and overlap
$\delta\approx d/4$. In LXZZ-setting runs, the same Helmholtz-harmonic basis
construction is instead placed on the LXZZ local partition so that both the
coarse basis and local solver see rectangular subdomains. The LXZZ local
partition uses `partitionMesh2D` on an aligned $H_{\rm LXZZ}^{-1}\times
H_{\rm LXZZ}^{-1}$ rectangular grid with $H_{\rm LXZZ}\le 1/\kappa$ and
`linearPartitionOfUnity2D` tensor-product linear weights, normalized nodally
inside the preconditioner builder. This replaces hat-support patches in the
Hu--Li economic cross-study because the economic trace parameterization
requires nondegenerate axis-aligned rectangular subdomain boundaries. Both

$$
\varepsilon=0
\quad\text{and}\quad
\varepsilon=\kappa
$$

are tested.

## Experiment Parameters

The fine-grid rule is

$$
h\approx\kappa^{-(2p+1)/(2p)}.
$$

The implementation uses unit asymptotic constants, then rounds upward so that
$h$, $d$, and $\delta$ align. Standard unrestarted GMRES uses zero initial
guess, relative Euclidean residual tolerance $10^{-6}$, and at most 1000
iterations.

Rows estimated above 200 GB require explicit permission. The gate uses sparse
matrix and sparse-LU storage models; it does not use dense $16r^2$ coarse
operator or full-LU storage estimates. Drivers checkpoint every case and report
whether it ran, was queued by the interactive runtime cap, requires permission,
or exceeds the hard memory limit.
