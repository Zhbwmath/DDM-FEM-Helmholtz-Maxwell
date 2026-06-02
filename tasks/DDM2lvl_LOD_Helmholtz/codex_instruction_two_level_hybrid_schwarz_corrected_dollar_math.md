# Codex Instruction: Corrected Reproduction Plan for
# “Two-level hybrid Schwarz Preconditioners for The Helmholtz Equation with high wave number”

Target article:

- **Title:** Two-level hybrid Schwarz Preconditioners for The Helmholtz Equation with high wave number
- **arXiv:** 2408.07669
- **Version used:** v2, revised 25 February 2025
- **Authors:** Peipei Lu, Xuejun Xu, Bowen Zheng, Jun Zou
- **URL:** https://arxiv.org/abs/2408.07669

This file corrects a crucial implementation issue:

> The paper defines $Q_m^{(1)}$ and $Q_m^{(2)}$ as **preconditioned operators on finite-element functions**, i.e. $Q_m^{(i)}=B_i^{-1}A$, not directly as residual-to-correction preconditioners $B_i^{-1}$.

Therefore, in left-preconditioned GMRES, distinguish


$$
M_0^{-1}
\quad\text{from}\quad
Q_0=M_0^{-1}A.
$$


Do **not** update a residual by

```text
e = r - Q0(r)
```

when `r` is a residual. The correct residual update is


$$
r_c=r-AM_0^{-1}r.
$$


---

## 1. Model Problem

Solve


$$
-\Delta u-\kappa^2u=f\quad\text{in }\Omega=(0,1)^2,
$$


with impedance boundary condition


$$
\partial_nu-i\kappa u=g\quad\text{on }\Gamma=\partial\Omega.
$$


Use the exact plane wave


$$
u(x,y)=\exp\!\left(i\kappa\frac{x+y}{\sqrt2}\right)
      =\exp(i\kappa d\cdot x),
\qquad
d=(1/\sqrt2,1/\sqrt2).
$$


Then


$$
f=0,
$$


and


$$
g=\partial_nu-i\kappa u
  =
  i\kappa(d\cdot n-1)u.
$$


---

## 2. Fine Finite-Element Discretization

Use conforming $P_1$ finite elements on a uniform triangular mesh of the unit square.

Let


$$
V_h=\operatorname{span}\{\phi_1,\ldots,\phi_n\}.
$$


For coefficient vectors $u,v\in\mathbb C^n$,


$$
u_h=\sum_j u_j\phi_j,\qquad v_h=\sum_j v_j\phi_j.
$$


Use the convention


$$
a(u_h,v_h)=v^{\mathrm H}Au.
$$


The Helmholtz form is


$$
a(u,v)
=
\int_\Omega \nabla u\cdot\nabla\bar v
-\kappa^2\int_\Omega u\bar v
-i\kappa\int_\Gamma u\bar v.
$$


Thus


$$
A=S-\kappa^2M-i\kappa N,
$$


where


$$
S_{ij}=\int_\Omega \nabla\phi_j\cdot\nabla\bar\phi_i,\qquad
M_{ij}=\int_\Omega \phi_j\bar\phi_i,\qquad
N_{ij}=\int_\Gamma \phi_j\bar\phi_i.
$$


The linear system is


$$
Au=F.
$$


The energy inner-product matrix is


$$
D_\kappa=S+\kappa^2M,
$$


so that


$$
(u_h,v_h)_{1,\kappa}=v^{\mathrm H}D_\kappa u.
$$


---

## 3. Mesh Sizes

The paper uses


$$
h\sim\kappa^{-3/2}.
$$


Implement this by

```text
n_h = rounding_rule(C_h * kappa^(3/2))
h = 1 / n_h
```

with configurable

```text
C_h
rounding_rule in {ceil, round, floor_to_valid_refinement}
```

Default:

```text
C_h = 1
rounding_rule = ceil
```

For experiments with prescribed $h=2^{-j}$, directly set the uniform refinement level.

---

## 4. Coarse Space and Quasi-Interpolation

Let


$$
V_H=\operatorname{span}\{\Phi_1,\ldots,\Phi_{N_H}\}\subset V_h.
$$


Let


$$
P\in\mathbb C^{n\times N_H}
$$


be the coarse-to-fine injection matrix. The $p$-th column of $P$ is the fine-grid coefficient vector of $\Phi_p$.

The quasi-interpolation is


$$
I_H v
=
\sum_{p=1}^{N_H}
\frac{(v,\Phi_p)_{L^2(\Omega)}}{(1,\Phi_p)_{L^2(\Omega)}}\Phi_p.
$$


Let


$$
G=\operatorname{diag}\left((1,\Phi_p)_{L^2(\Omega)}\right)_{p=1}^{N_H}.
$$


For a fine-grid coefficient vector $v$, the coarse coefficient vector of $I_Hv_h$ is


$$
Jv=G^{-1}P^{\mathrm H}Mv.
$$


Thus


$$
W_h=\ker I_H=\ker J.
$$


---

## 5. Localized LOD Correctors

For every coarse element $T$, define the $m$-layer patch $\omega^m(T)$ and


$$
W_{h,T,m}
=
\{w_h\in W_h:\operatorname{supp}w_h\subset \omega^m(T)\}.
$$


Let


$$
Z_T\in\mathbb C^{n\times n_T}
$$


be a basis matrix for $W_{h,T,m}$.

Let $A_T$ be the matrix of the element-local form $a_T(\cdot,\cdot)$, extended by zero to the global fine space. For the volume Helmholtz term,


$$
a_T(u,v)=
\int_T\nabla u\cdot\nabla\bar v
-\kappa^2\int_Tu\bar v.
$$


Check boundary-element conventions against the paper if adding boundary terms locally.

### 5.1 Primal Corrector

For $v_H=P\alpha$, define $C_{T,m}v_H\in W_{h,T,m}$ by


$$
a(C_{T,m}v_H,w)=a_T(v_H,w)
\qquad\forall w\in W_{h,T,m}.
$$


Write $C_{T,m}v_H=Z_T\gamma$. Then


$$
Z_T^{\mathrm H}AZ_T\gamma
=
Z_T^{\mathrm H}A_TP\alpha.
$$


Hence


$$
C_{T,m}
=
Z_T
\left(Z_T^{\mathrm H}AZ_T\right)^{-1}
Z_T^{\mathrm H}A_TP.
$$


Here $C_{T,m}\in\mathbb C^{n\times N_H}$.

### 5.2 Adjoint Corrector

The adjoint corrector satisfies


$$
a(w,C_{T,m}^{\ast}v_H)=a_T(w,v_H)
\qquad\forall w\in W_{h,T,m}.
$$


Write $C_{T,m}^{\ast}v_H=Z_T\gamma^\ast$. Then


$$
Z_T^{\mathrm H}A^{\mathrm H}Z_T\gamma^\ast
=
Z_T^{\mathrm H}A_T^{\mathrm H}P\alpha.
$$


Hence


$$
C_{T,m}^{\ast}
=
Z_T
\left(Z_T^{\mathrm H}A^{\mathrm H}Z_T\right)^{-1}
Z_T^{\mathrm H}A_T^{\mathrm H}P.
$$


### 5.3 Global Corrector Matrices

Set


$$
C_m=\sum_T C_{T,m},
\qquad
C_m^\ast=\sum_T C_{T,m}^\ast.
$$


Both are $n\times N_H$ matrices.

---

## 6. LOD Trial and Test Basis Matrices

The localized trial and test spaces are


$$
V_{H,m}=(I-C_m)V_H,\qquad
V_{H,m}^\ast=(I-C_m^\ast)V_H.
$$


Since $C_m$ and $C_m^\ast$ map coarse coefficients into fine-grid correctors, the practical matrices are


$$
B=P-C_m,
\qquad
B_\ast=P-C_m^\ast.
$$


Use

```text
B_trial = P - C_m
B_test  = P - C_m_star
```

Do **not** use

```text
B_trial = (I - C_m) * P
```

unless $C_m$ has been explicitly assembled as an $n\times n$ fine-space operator.

---

## 7. LOD Coarse Operator $Q_0$

The paper defines $Q_{0,m}:V_h\to V_{H,m}$ by


$$
a(Q_{0,m}v,w)=a(v,w)
\qquad\forall w\in V_{H,m}^{\ast}.
$$


For a fine vector $x$, seek


$$
Q_{0,m}x=Bc.
$$


Testing with $B_\ast y$,


$$
a(Bc,B_\ast y)=a(x,B_\ast y)
\qquad\forall y.
$$


Using $a(u,v)=v^{\mathrm H}Au$,


$$
B_\ast^{\mathrm H}ABc=B_\ast^{\mathrm H}Ax.
$$


Define


$$
A_H^{\rm LOD}=B_\ast^{\mathrm H}AB.
$$


Then


$$
c=(A_H^{\rm LOD})^{-1}B_\ast^{\mathrm H}Ax.
$$


Therefore,


$$
Q_0x
=
B(A_H^{\rm LOD})^{-1}B_\ast^{\mathrm H}Ax.
$$


Define the coarse residual-to-correction inverse


$$
M_0^{-1}
=
B(A_H^{\rm LOD})^{-1}B_\ast^{\mathrm H}.
$$


Then


$$
Q_0=M_0^{-1}A.
$$


Implementation:

```text
function z0 = apply_M0_inverse(r):
    rhs_H = B_test^H * r
    coeff = A_H_LOD \ rhs_H
    z0 = B_trial * coeff
    return z0

function q0 = apply_Q0(x):
    return apply_M0_inverse(A * x)
```

---

## 8. Paper’s Adjoint Coarse Operator $Q_0^\ast$

The paper also defines $Q_{0,m}^\ast:V_h\to V_{H,m}^\ast$ by


$$
a(w,Q_{0,m}^{\ast}v)=a(w,v)
\qquad\forall w\in V_{H,m}.
$$


Let


$$
Q_{0,m}^{\ast}x=B_\ast d.
$$


Testing with $By$,


$$
B^{\mathrm H}A^{\mathrm H}B_\ast d
=
B^{\mathrm H}A^{\mathrm H}x.
$$


Since


$$
B^{\mathrm H}A^{\mathrm H}B_\ast=(A_H^{\rm LOD})^{\mathrm H},
$$


we get


$$
d=(A_H^{\rm LOD})^{-\mathrm H}B^{\mathrm H}A^{\mathrm H}x.
$$


Thus


$$
Q_0^\ast x
=
B_\ast(A_H^{\rm LOD})^{-\mathrm H}B^{\mathrm H}A^{\mathrm H}x.
$$


Important: this $Q_0^\ast$ is **not** the same as the hybrid transpose $(I-Q_0)^T$.

---

## 9. Energy Adjoint in the Hybrid Formula

The hybrid transpose $(I-Q_0)^T$ is the adjoint with respect to $(\cdot,\cdot)_{1,\kappa}$.

For any matrix $T$,


$$
T^{T_D}=D_\kappa^{-1}T^{\mathrm H}D_\kappa.
$$


Hence


$$
(I-Q_0)^T
=
(I-Q_0)^{T_D}
=
D_\kappa^{-1}(I-Q_0)^{\mathrm H}D_\kappa.
$$


The Euclidean adjoint of $Q_0$ is


$$
Q_0^{\mathrm H}
=
A^{\mathrm H}B_\ast(A_H^{\rm LOD})^{-\mathrm H}B^{\mathrm H}.
$$


Implementation:

```text
function y = Q0_Euclidean_adjoint(w):
    rhs = B_trial^H * w
    coeff = A_H_LOD^H \ rhs
    y = A^H * B_test * coeff
    return y

function y = energy_adjoint_I_minus_Q0(z):
    w = D_kappa * z
    t = w - Q0_Euclidean_adjoint(w)
    y = solve(D_kappa, t)
    return y
```

---

## 10. Dirichlet-Type Local Solver

For the first preconditioner, the local space is


$$
\widetilde V_{h,\ell}
=
\{v_h|_{\Omega_\ell}:v_h=0\text{ on }\partial\Omega_\ell\setminus\Gamma\}.
$$


The paper defines


$$
a_\ell(Q_\ell v_h,w_{h,\ell})
=
a(v_h,w_{h,\ell})
\qquad
\forall w_{h,\ell}\in\widetilde V_{h,\ell}.
$$


Let $E_\ell^D$ extend local Dirichlet degrees of freedom into global degrees of freedom. Let


$$
A_{\ell,D}
=
(E_\ell^D)^{\mathrm H}A_\ell E_\ell^D.
$$


Define the local residual inverse


$$
M_{\ell,D}^{-1}
=
E_\ell^D A_{\ell,D}^{-1}(E_\ell^D)^{\mathrm H}.
$$


Then


$$
Q_\ell=M_{\ell,D}^{-1}A.
$$


Set


$$
M_D^{-1}=\sum_{\ell=1}^N M_{\ell,D}^{-1}.
$$


Implementation:

```text
function z = apply_MD_inverse(r):
    z = 0
    for ell in subdomains:
        rloc = E_D_ell^H * r
        yloc = A_D_ell \ rloc
        z += E_D_ell * yloc
    return z
```

---

## 11. First Hybrid Operator $Q_m^{(1)}$

The paper defines


$$
Q_m^{(1)}
=
Q_0
+
(I-Q_0)^T
\sum_{\ell=1}^N Q_\ell
(I-Q_0).
$$


Matrix form:


$$
Q^{(1)}
=
Q_0
+
(I-Q_0)^{T_D}
M_D^{-1}A
(I-Q_0).
$$


This acts on function vectors $x$. Function-level implementation:

```text
function y = apply_Q1_operator(x):
    q0 = apply_Q0(x)
    e  = x - q0
    s  = apply_MD_inverse(A * e)
    y  = q0 + energy_adjoint_I_minus_Q0(s)
    return y
```

---

## 12. Residual-Level Preconditioner $B_1^{-1}$

For actual left-preconditioned GMRES, implement


$$
B_1^{-1}
=
M_0^{-1}
+
(I-Q_0)^{T_D}
M_D^{-1}
(I-AM_0^{-1}).
$$


For a residual $r$,


$$
z_0=M_0^{-1}r,
$$


$$
r_c=r-Az_0,
$$


$$
s=M_D^{-1}r_c,
$$


$$
z=z_0+(I-Q_0)^{T_D}s.
$$


Implementation:

```text
function z = apply_B1_inverse(r):
    z0 = apply_M0_inverse(r)
    rc = r - A * z0
    s  = apply_MD_inverse(rc)
    z  = z0 + energy_adjoint_I_minus_Q0(s)
    return z
```

Then


$$
B_1^{-1}A=Q_m^{(1)}.
$$


---

## 13. Impedance-Type Local Solver

The second local solver uses


$$
c_\ell(u,v)
=
\int_{\Omega_\ell}\nabla u\cdot\nabla\bar v
+
\kappa^2\int_{\Omega_\ell}u\bar v
-
i\kappa\int_{\partial\Omega_\ell}u\bar v.
$$


The sign is important: the local form has $+\kappa^2$, not $-\kappa^2$.

The paper defines $P_\ell$ by


$$
c_\ell(P_\ell v,w_{h,\ell})
=
a(v,\Pi_h\chi_\ell w_{h,\ell})
\qquad
\forall w_{h,\ell}\in V_{h,\ell}.
$$


Let $E_\ell$ extend all local subdomain degrees of freedom to global degrees of freedom. Let $X_\ell$ be the diagonal matrix of nodal values of the partition-of-unity function $\chi_\ell$. Define


$$
W_\ell=E_\ell X_\ell.
$$


Then $W_\ell y$ is the global coefficient vector of $\Pi_h(\chi_\ell w_{h,\ell})$.

Let $C_\ell$ be the local matrix for $c_\ell$. For a function vector $x$,


$$
C_\ell p_\ell=W_\ell^{\mathrm H}Ax.
$$


Hence


$$
\Pi_h\chi_\ell P_\ell x
=
W_\ell C_\ell^{-1}W_\ell^{\mathrm H}Ax.
$$


Define


$$
M_{\ell,I}^{-1}
=
W_\ell C_\ell^{-1}W_\ell^{\mathrm H}.
$$


Set


$$
M_I^{-1}=\sum_{\ell=1}^N M_{\ell,I}^{-1}.
$$


Implementation:

```text
function z = apply_MI_inverse(r):
    z = 0
    for ell in subdomains:
        rloc = W_ell^H * r
        yloc = C_ell \ rloc
        z += W_ell * yloc
    return z
```

---

## 14. Second Hybrid Operator $Q_m^{(2)}$

The paper defines


$$
Q_m^{(2)}
=
Q_0
+
(I-Q_0)^T
\sum_{\ell=1}^N
\Pi_h\chi_\ell P_\ell
(I-Q_0).
$$


Matrix form:


$$
Q^{(2)}
=
Q_0
+
(I-Q_0)^{T_D}
M_I^{-1}A
(I-Q_0).
$$


Function-level implementation:

```text
function y = apply_Q2_operator(x):
    q0 = apply_Q0(x)
    e  = x - q0
    s  = apply_MI_inverse(A * e)
    y  = q0 + energy_adjoint_I_minus_Q0(s)
    return y
```

---

## 15. Residual-Level Preconditioner $B_2^{-1}$

For left-preconditioned GMRES, implement


$$
B_2^{-1}
=
M_0^{-1}
+
(I-Q_0)^{T_D}
M_I^{-1}
(I-AM_0^{-1}).
$$


For a residual $r$,


$$
z_0=M_0^{-1}r,
$$


$$
r_c=r-Az_0,
$$


$$
s=M_I^{-1}r_c,
$$


$$
z=z_0+(I-Q_0)^{T_D}s.
$$


Implementation:

```text
function z = apply_B2_inverse(r):
    z0 = apply_M0_inverse(r)
    rc = r - A * z0
    s  = apply_MI_inverse(rc)
    z  = z0 + energy_adjoint_I_minus_Q0(s)
    return z
```

Then


$$
B_2^{-1}A=Q_m^{(2)}.
$$


---

## 16. Standard $P_1$ Coarse Solver for Comparison

For Tables 5.7 and 5.8, replace the LOD coarse solver by the standard $P_1$-FEM coarse solver.

Define


$$
A_H=P^{\mathrm H}AP,
$$


$$
M_H^{-1}=P A_H^{-1}P^{\mathrm H},
$$


$$
Q_H=M_H^{-1}A.
$$


The comparison residual-level preconditioners are


$$
\widetilde B_1^{-1}
=
M_H^{-1}
+
(I-Q_H)^{T_D}
M_D^{-1}
(I-AM_H^{-1}),
$$


and


$$
\widetilde B_2^{-1}
=
M_H^{-1}
+
(I-Q_H)^{T_D}
M_I^{-1}
(I-AM_H^{-1}).
$$


---

## 17. GMRES Settings

The paper uses standard GMRES with Euclidean residual minimization.

Use:

```text
initial_guess = zero
relative_residual_tolerance = 1e-6
```

Prefer unrestarted GMRES. If restarted GMRES is necessary, record the restart value.

Recommended solver:

```text
GMRES(A, F, left_preconditioner = apply_Bi_inverse)
```

where `apply_Bi_inverse` is one of:

```text
apply_B1_inverse
apply_B2_inverse
apply_B1tilde_inverse
apply_B2tilde_inverse
```

Record


$$
\frac{\|F-Au_j\|_2}{\|F-Au_0\|_2}.
$$


---

## 18. Experiment 1: Dependence on Wave Number

### 18.1 $Q_m^{(1)}$

Use


$$
h\sim\kappa^{-3/2},
\qquad
H=\delta=H_{\rm sub}/2\sim\kappa^{-1},
\qquad
m=\log_2(\kappa)-1.
$$


Run

```text
kappa_values = [16, 32, 64, 128, 256, 500]
```

Target table:

| $\kappa$ | Iterations |
|---:|---:|
| 16  | 9 |
| 32  | 8 |
| 64  | 8 |
| 128 | 8 |
| 256 | 8 |
| 500 | 8 |

### 18.2 $Q_m^{(2)}$

Use


$$
h\sim\kappa^{-3/2},
\qquad
H\sim\kappa^{-1},
\qquad
H_{\rm sub}=4H,
\qquad
\delta=2H,
\qquad
m=\log_2(\kappa)-1.
$$


Run the same `kappa_values`.

Target table:

| $\kappa$ | Iterations |
|---:|---:|
| 16  | 7 |
| 32  | 7 |
| 64  | 7 |
| 128 | 7 |
| 256 | 7 |
| 500 | 7 |

---

## 19. Experiment 2: $h$-Independence

Fix

```text
kappa = 80
H ~ kappa^(-1)
m = 2
```

### 19.1 $Q_m^{(1)}$

Use


$$
H=\delta=H_{\rm sub}/2.
$$


Run

```text
h_values = [2^(-10), 2^(-11), 2^(-12), 2^(-13)]
```

Target table:

| $h$ | Iterations |
|---:|---:|
| $2^{-10}$ | 10 |
| $2^{-11}$ | 9 |
| $2^{-12}$ | 9 |
| $2^{-13}$ | 9 |

### 19.2 $Q_m^{(2)}$

Use


$$
\delta=2H,
\qquad
H_{\rm sub}=4H.
$$


Run the same `h_values`.

Target table:

| $h$ | Iterations |
|---:|---:|
| $2^{-10}$ | 9 |
| $2^{-11}$ | 8 |
| $2^{-12}$ | 8 |
| $2^{-13}$ | 8 |

---

## 20. Experiment 3: Oversampling Parameter $m$

Fix

```text
kappa = 128
```

### 20.1 $Q_m^{(1)}$

Use


$$
h\sim\kappa^{-3/2},
\qquad
H=\delta=H_{\rm sub}/2=\kappa^{-1}.
$$


Run

```text
m_values = [6, 5, 4, 3, 2, 1]
```

Target table:

| $m$ | Iterations |
|---:|---:|
| 6 | 8 |
| 5 | 8 |
| 4 | 8 |
| 3 | 8 |
| 2 | 9 |
| 1 | 11 |

### 20.2 $Q_m^{(2)}$

Use


$$
h\sim\kappa^{-3/2},
\qquad
H=\kappa^{-1},
\qquad
H_{\rm sub}=4H,
\qquad
\delta=2H.
$$


Run the same `m_values`.

Target table:

| $m$ | Iterations |
|---:|---:|
| 6 | 7 |
| 5 | 7 |
| 4 | 7 |
| 3 | 7 |
| 2 | 8 |
| 1 | 10 |

---

## 21. Experiment 4: Subdomain Size and Overlap

Let $H_0=H_0(\kappa)$ satisfy


$$
H_0\kappa\sim1.
$$


Use


$$
h\sim\kappa^{-3/2},
\qquad
m=2.
$$


Also implement $\widetilde Q^{(1)}$ and $\widetilde Q^{(2)}$ using the standard $P_1$-FEM coarse solver.

Run

```text
kappa_values = [40, 80, 120, 160]
```

### 21.1 $Q_m^{(1)}$ and $\widetilde Q^{(1)}$

Use


$$
H=\delta=H_{\rm sub}/2.
$$


Target table. Parentheses are standard $P_1$-coarse results.

| $\kappa$ | $H_{\rm sub}=2H_0$ | $H_{\rm sub}=H_0$ |
|---:|---:|---:|
| 40  | 9 (25)     | 8 (23) |
| 80  | 9 (55)     | 8 (47) |
| 120 | 9 (>100)   | 8 (85) |
| 160 | 9 (>100)   | 8 (>100) |

### 21.2 $Q_m^{(2)}$ and $\widetilde Q^{(2)}$

Use


$$
H=H_0,
\qquad
\delta=H_{\rm sub}/2.
$$


Target table. Parentheses are standard $P_1$-coarse results.

| $\kappa$ | $\delta=2H_0$ | $\delta=4H_0$ |
|---:|---:|---:|
| 40  | 7 (26)     | 6 (21) |
| 80  | 7 (49)     | 6 (43) |
| 120 | 7 (89)     | 6 (77) |
| 160 | 6 (>100)   | 6 (>100) |

---

## 22. Experiment 5: Small-Overlap Behavior for $Q_m^{(1)}$

Use


$$
h\sim\kappa_{\max}^{-3/2},
\qquad
\kappa_{\max}=120,
$$


and keep this fine mesh fixed for all $\kappa$.

Choose


$$
H=H_0=H_{\rm sub}/2,
\qquad
H_0\kappa\sim1,
\qquad
m=\log_2(\kappa).
$$


Run

```text
kappa_values = [40, 80, 120]
delta_values = [H0, 4h, 2h, h]
```

Target table:

| $\kappa$ | $\delta=H_0$ | $\delta=4h$ | $\delta=2h$ | $\delta=h$ |
|---:|---:|---:|---:|---:|
| 40  | 9 | 11 | 12 | 13 |
| 80  | 8 | 12 | 14 | 15 |
| 120 | 8 | 19 | 22 | 26 |

---

## 23. Required Outputs

Generate:

1. A script or package that runs all experiments.
2. A table generator for Tables 5.1–5.9.
3. A log file for each run containing:
   - $\kappa$,
   - $h$,
   - $H$,
   - $H_{\rm sub}$,
   - $\delta$,
   - $m$,
   - global degrees of freedom,
   - coarse degrees of freedom,
   - number of subdomains,
   - local degrees-of-freedom statistics,
   - GMRES iterations,
   - final relative residual,
   - LOD corrector construction time,
   - coarse solve time,
   - local solver setup time,
   - local preconditioner application time,
   - total runtime.
4. A reproducibility note explaining:
   - constants used in $h\sim\kappa^{-3/2}$,
   - constants used in $H\sim\kappa^{-1}$,
   - rounding rules,
   - subdomain and overlap construction,
   - implementation of the $D_\kappa$-adjoint,
   - whether GMRES applies $B_i^{-1}$ as a left preconditioner or explicitly applies $Q_m^{(i)}$.
5. A comparison report showing deviations from the target iteration counts.

---

## 24. Debugging Checklist

### 24.1 LOD Correctors

Verify:


$$
J C_{T,m}\approx0,
\qquad
J C_{T,m}^\ast\approx0.
$$


Verify:


$$
B=P-C_m,
\qquad
B_\ast=P-C_m^\ast.
$$


Verify:


$$
A_H^{\rm LOD}=B_\ast^{\mathrm H}AB.
$$


### 24.2 Coarse Solver

Verify:


$$
M_0^{-1}=B(A_H^{\rm LOD})^{-1}B_\ast^{\mathrm H},
\qquad
Q_0=M_0^{-1}A.
$$


### 24.3 Residual Update

For residual-level preconditioning, verify:


$$
r_c=r-AM_0^{-1}r.
$$


Do not use $r-Q_0r$ for residuals.

### 24.4 Energy Adjoint

Verify numerically, for random $x,y$,


$$
((I-Q_0)x,y)_{1,\kappa}
=
(x,(I-Q_0)^{T_D}y)_{1,\kappa}.
$$


### 24.5 Identity Test

For a small problem, explicitly assemble $B_i^{-1}$ and $Q_m^{(i)}$ and verify


$$
B_i^{-1}A\approx Q_m^{(i)}.
$$


### 24.6 Local Solvers

For $Q_m^{(1)}$, enforce homogeneous Dirichlet conditions on artificial subdomain boundaries.

For $Q_m^{(2)}$, verify the local matrix uses


$$
+\kappa^2M_\ell
$$


inside $c_\ell$, not $-\kappa^2M_\ell$.

### 24.7 GMRES

Verify:

```text
initial_guess = zero
relative_residual_tolerance = 1e-6
residual_norm = euclidean
```

Record whether GMRES is restarted.

---

## 25. Suggested Code Structure

```text
project/
  README.md
  config.yaml
  run_all_experiments.py
  src/
    mesh.py
    fem_assembly.py
    helmholtz_problem.py
    quasi_interpolation.py
    lod_correctors.py
    lod_coarse_solver.py
    energy_adjoint.py
    subdomains.py
    local_dirichlet_solver.py
    local_impedance_solver.py
    hybrid_preconditioners.py
    gmres_driver.py
    tables.py
    logging_utils.py
  experiments/
    exp1_wave_number.py
    exp2_h_independence.py
    exp3_oversampling.py
    exp4_subdomain_overlap.py
    exp5_small_overlap.py
  results/
    logs/
    tables/
    figures/
```

Suggested configuration:

```yaml
fine_mesh:
  C_h: 1.0
  rounding_rule: ceil

coarse_mesh:
  C_H: 1.0
  rounding_rule: ceil

lod:
  enforce_kernel_exactly: true
  use_adjoint_correctors: true

hybrid:
  residual_level_left_preconditioner: true
  use_D_kappa_adjoint: true

gmres:
  tolerance: 1.0e-6
  restart: null
  initial_guess: zero
  residual_norm: euclidean

output:
  save_logs: true
  save_tables: true
  save_timings: true
```
