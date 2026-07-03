# Codex Instruction: Three-Level LOD-DDM Preconditioner and Power-Contraction Diagnostics

## 0. Purpose

Implement and test a three-level Helmholtz preconditioner in which the exact global LOD coarse solve is replaced by a one-level additive Schwarz/DDM solver on the LOD coarse system.

The main goal is not only to record GMRES iteration counts, but also to verify whether the observed convergence is explained by a power-contraction mechanism:

$E_0 := I - M_0^{-1}A_0$,

$|E_0^s|_{G_0} < 1$

for some $s$ proportional to the effective coarse-subdomain graph diameter.

Here $A_0$ is the global LOD coarse matrix, $M_0^{-1}$ is the one-level additive coarse solver, and $G_0$ is a positive definite discrete $H^1_k$-type energy matrix on the LOD coarse space.

The implementation should be modular. First implement the coarse-level diagnostics. Then integrate the coarse solver into the full hybrid Schwarz/LOD preconditioner.

---

## 1. Continuous and Discrete Problem

Consider the Helmholtz-type variational form

### $a_\varepsilon(u,v) = (\nabla u,\nabla v)_\Omega (k^2+i\varepsilon)(u,v)_\Omega i\sqrt{k^2+i\varepsilon}(u,v)_{\Gamma_R}$,

where:

* $k$ is the wave number;
* $\varepsilon \ge 0$ is the artificial absorption parameter;
* $\Gamma_R$ is the impedance boundary part;
* the square root is taken with positive real part.

The non-absorbing case is $\varepsilon=0$. Some experiments should also use $\varepsilon=k$.

Let $V_h$ be the fine finite element space and $V_H$ the coarse nodal finite element space. Let

$V_H = \operatorname{span}{\Phi_i : i\in I_H}$.

Let $I_H:V_h\to V_H$ be a projective quasi-interpolation operator, and define the fine-scale kernel

$W_h := \ker I_H$.

The localized LOD correction $C_m:V_H\to W_h$ is defined patchwise. The corrected global LOD trial basis is

$\psi_i := (I-C_m)\Phi_i$.

If an adjoint/Petrov LOD formulation is used, also build the adjoint corrected test basis

$\psi_i^\star := (I-C_m^\star)\Phi_i$.

The global LOD coarse space is

$V_{H,m}^{\mathrm{LOD}} := \operatorname{span}{\psi_i}$.

The test space is

$V_{H,m}^{\star,\mathrm{LOD}} := \operatorname{span}{\psi_i^\star}$.

---

## 2. Coarse Subdomains and Local LOD Subproblems

Let ${\Omega_{0,\ell}}_{\ell=1}^{N_0}$ be the coarse-level subdomain partition. Let $H_{0,\ell}$ denote the diameter of $\Omega_{0,\ell}$. In the target regime, $H_{0,\ell}$ should be independent of $k$.

For each $\Omega_{0,\ell}$, define an oversampled LOD domain $\widetilde\Omega_{0,\ell}$ by collecting the supports of corrected LOD basis functions whose original coarse support belongs to $\overline{\Omega}_{0,\ell}$:

### $\widetilde\Omega_{0,\ell} \bigcup{\operatorname{supp}(C_m\Phi_i): \operatorname{supp}\Phi_i\subset \overline{\Omega}_{0,\ell}}$.

Let $\chi_{0,\ell}$ be a coarse-level partition of unity satisfying

$\sum_\ell \chi_{0,\ell}=1$,

$0\le \chi_{0,\ell}\le 1$,

$|\nabla \chi_{0,\ell}|_{L^\infty}\lesssim \delta_{0,\ell}^{-1}$.

Let $\chi_{0,\ell}^{>}$ be an enlarged cutoff satisfying

$0\le \chi_{0,\ell}^{>}\le 1$,

$\operatorname{supp}\chi_{0,\ell}^{>}\subset \Omega_{0,\ell}$,

$\chi_{0,\ell}^{>}=1$ on $\operatorname{supp}\chi_{0,\ell}$,

$|\nabla\chi_{0,\ell}^{>}|_{L^\infty}\lesssim \delta_{0,\ell}^{-1}$.

The local fine space is

$V_{h,\ell}:=V_h(\widetilde\Omega_{0,\ell})$.

The local coarse space is

$V_{H,\ell}:=V_H(\widetilde\Omega_{0,\ell})$.

Let $I_{H,\ell}:V_{h,\ell}\to V_{H,\ell}$ be the local quasi-interpolation, and define

$W_{h,\ell}:=\ker I_{H,\ell}$.

The local LOD correction $C_{\ell,m}$ is the localized analogue of $C_m$. The local corrected space is

$V_{H,m,\ell}:=(I-C_{\ell,m})V_{H,\ell}$.

When the local and global correctors agree on interior patches, use the global corrected basis functions $\psi_i=(I-C_m)\Phi_i$ to assemble the local systems.

---

## 3. Global LOD Coarse Matrix

Let $\Psi$ be the matrix whose columns are the global corrected trial basis functions $\psi_i$ in fine-grid coordinates.

Let $\Psi^\star$ be the matrix whose columns are the global corrected adjoint test functions $\psi_i^\star$.

If a symmetric/Galerkin version is used, set $\Psi^\star=\Psi$.

Let $A_h$ be the fine-grid Helmholtz matrix corresponding to $a_\varepsilon(\cdot,\cdot)$.

Define the global LOD coarse matrix

$A_0 := (\Psi^\star)^T A_h \Psi$.

Equivalently,

$(A_0)_{ij}=a_\varepsilon(\psi_j,\psi_i^\star)$.

The exact global LOD coarse solve is

$A_0 z = r_0$.

This exact solve is the baseline coarse solve used in the two-level method. The three-level method replaces this exact solve by a one-level additive Schwarz/DDM solver on the coarse matrix $A_0$.

---

## 4. Local Coarse Solver

For each coarse subdomain $\Omega_{0,\ell}$, define a local corrected basis matrix $\Psi_\ell$ and a local adjoint corrected test basis matrix $\Psi_\ell^\star$.

The local coarse matrix is

$A_{0,\ell}:=(\Psi_\ell^\star)^T A_{h,\ell}\Psi_\ell$,

where $A_{h,\ell}$ is the local Helmholtz matrix on $\widetilde\Omega_{0,\ell}$, including the local impedance boundary term on artificial boundaries.

The variational definition of the local coarse solver is:

Find $Q_{0,m,\ell}v_h\in V_{H,m,\ell}$ such that

### $a_\ell(Q_{0,m,\ell}v_h,v_{H,m,\ell}^\star) = a(v_h,(I-C_m^\star)\Pi_H\chi_{0,\ell}^{>}I_{H,\ell}v_{H,m,\ell}^\star)$

for all $v_{H,m,\ell}^\star\in V_{H,m,\ell}^{\star}$.

In algebraic form, this means:

1. Form the local right-hand side using the enlarged cutoff $\chi_{0,\ell}^{>}$.
2. Solve the local LOD coarse system with $A_{0,\ell}$.
3. Extend/reconstruct the local solution into the global LOD coarse coefficient vector using the partition-of-unity weight $\chi_{0,\ell}$.

Implement this as a function:

`apply_M0inv(r0)`

where `r0` is a global LOD coarse residual vector.

The additive one-level coarse solver is

### $M_0^{-1}r_0 = \sum_{\ell=1}^{N_0} D_\ell R_\ell^T A_{0,\ell}^{-1} R_\ell^{>} r_0$.

Here:

* $R_\ell^{>}$ is the enlarged-cutoff restriction induced by $\chi_{0,\ell}^{>}$;
* $R_\ell^T$ is local-to-global extension;
* $D_\ell$ is the partition-of-unity reconstruction induced by $\chi_{0,\ell}$;
* $A_{0,\ell}^{-1}$ should be applied by a direct sparse solver in the first implementation.

Important: `R_l_greater` and `D_l` should not be treated as ordinary Boolean restrictions unless that is exactly how the finite element code realizes the cutoff functions. Prefer implementing them through the finite element matrices induced by $\chi_{0,\ell}^{>}$ and $\chi_{0,\ell}$.

---

## 5. Coarse Error Propagation Operator

After implementing `apply_M0inv`, define the coarse preconditioned operator

$S_0:=M_0^{-1}A_0$,

and the coarse error-propagation operator

$E_0:=I-S_0$.

For small coarse systems, explicitly assemble

$E_0 = I - M_0^{-1}A_0$.

For large systems, implement function handles:

`apply_S0(x) = apply_M0inv(A0*x)`

`apply_E0(x) = x - apply_S0(x)`.

---

## 6. Positive Definite Coarse Energy Norm

Do not use the indefinite Helmholtz matrix as the norm matrix.

Build a positive definite coarse energy matrix

$G_0 := \Psi^*(K_h+k^2M_h)\Psi$,

where:

* $K_h$ is the fine stiffness matrix;
* $M_h$ is the fine mass matrix;
* $\Psi$ is the global corrected LOD trial basis matrix.

Thus

$|x|_{G_0}^2 := x^*G_0x$

is the discrete LOD version of

### $|v|_{1,k}^2 = |\nabla v|_{L^2(\Omega)}^2 + k^2|v|_{L^2(\Omega)}^2$.

If $G_0$ is numerically singular, remove dependent corrected basis functions or use a stable QR/SVD compression of $\Psi$.

---

## 7. Main Numerical Diagnostics

### Diagnostic A: coarse power norm

For $s=1,\dots,s_{\max}$, compute

### $|E_0^s|_{G_0}=\sqrt{\lambda_{\max}\left((E_0^s)^*G_0E_0^s,,G_0\right)}$.

Implementation:

1. Form $E_0^s$ if feasible.
2. Solve the generalized eigenvalue problem

$(E_0^s)^*G_0E_0^s x = \lambda G_0x$.

3. Record

`norm_E_power[s] = sqrt(real(lambda_max))`.

The key quantity is

$s_0 := \min{s:|E_0^s|_{G_0}<1}$.

Check whether $s_0$ scales with the effective coarse-subdomain graph diameter rather than with $k$.

---

### Diagnostic B: field of values of the polynomially transformed operator

Define

$P_s:=I-E_0^s$.

Compute

### $\alpha_s=\lambda_{\min}\left(\frac{G_0P_s+P_s^*G_0}{2},G_0\right)$.

If $\alpha_s>0$, then $P_s$ has a positive field-of-values lower bound in the $G_0$ inner product.

Record:

* `alpha_s`;
* the first $s$ for which `alpha_s > 0`;
* comparison between this $s$ and the first $s$ for which $|E_0^s|_{G_0}<1$.

Expected relation:

If $|E_0^s|_{G_0}\le q<1$, then $\alpha_s\ge 1-q$.

---

### Diagnostic C: GMRES for the coarse preconditioned system

Run GMRES on

$S_0 z = M_0^{-1}r_0$

or equivalently use left-preconditioned GMRES for

$A_0 z=r_0$

with preconditioner $M_0^{-1}$.

Record:

* iteration count;
* residual history;
* convergence tolerance;
* $k$;
* $H$;
* $H_0$;
* $m$;
* $\varepsilon$;
* overlap width;
* number of coarse subdomains;
* effective graph diameter.

Compare observed GMRES convergence with the bound predicted by $|E_0^s|_{G_0}$.

---

### Diagnostic D: full three-level GMRES

Use the existing two-level LOD hybrid Schwarz preconditioner code as the baseline.

Do not rewrite the local fine-level solver unless necessary.

Replace only the exact global LOD coarse solve

`z = A0 \ r0`

by

`z = apply_M0inv(r0)`.

This gives the three-level preconditioner.

Run GMRES for the original fine-grid system

$A_h u_h=f_h$

with:

1. exact two-level preconditioner;
2. three-level preconditioner with one-level coarse solver;
3. optionally, three-level preconditioner with $s$ inner coarse Schwarz sweeps.

Record iteration counts and residual histories.

---

## 8. Optional: $s$-Sweep Coarse Solver

Define the $s$-sweep approximate inverse

### $G_0^{(s)}=\sum_{j=0}^{s-1}E_0^jM_0^{-1}$.

It satisfies the exact identity

$I-G_0^{(s)}A_0=E_0^s$.

Implement:

```text
function y = apply_G0_s(r0, s)
    y = 0
    e = r0
    for j = 0:s-1
        delta = apply_M0inv(e)
        y = y + delta
        e = e - A0*delta
    end
end
```

Use this only as a diagnostic or theoretical bridge. The main practical method may still use one coarse sweep inside outer GMRES.

Test whether the hybrid preconditioner using `apply_G0_s` recovers the behavior of the exact two-level method as $s$ increases.

---

## 9. Optional: LOD-Projected Trace-Transfer Diagnostic

This diagnostic is more advanced. Implement only after Diagnostics A-D work.

### 9.1 Directed interface set

Create the set of directed coarse interfaces

$\mathcal I={e=(\ell\to j):\Omega_{0,\ell}\text{ and }\Omega_{0,j}\text{ are neighbors}}$.

For each directed interface $e=(\ell\to j)$, define the impedance trace

$\gamma_e u = \partial_{n_\ell}u-i\kappa u$ on $\Gamma_{\ell j}$.

### 9.2 LOD trace matrix

Build a matrix $Z$ mapping LOD coarse coefficients to discrete interface traces:

$u_0 \mapsto Zu_0$.

Each column of $Z$ is the interface trace of one corrected LOD basis function.

If impedance flux extraction is difficult, start with the Dirichlet trace matrix as a preliminary diagnostic, then add the impedance trace.

### 9.3 Trace norm matrix

Build a block diagonal trace norm matrix $N_\Lambda$.

First implementation:

$N_\Lambda =$ interface mass matrix.

More refined implementation:

Use the local LOD-harmonic extension norm.

### 9.4 Projection onto LOD-generated trace space

Define

### $P_\Xi=Z(Z^*N_\Lambda Z)^\dagger Z^*N_\Lambda$.

Use a pseudoinverse or remove dependent trace columns if $Z^*N_\Lambda Z$ is rank-deficient.

### 9.5 Transfer matrix

If feasible, build the full interface transfer matrix $T$ by local homogeneous Helmholtz solves:

incoming impedance data on one interface $\to$ outgoing impedance data on neighboring interfaces.

Then compare:

$||T^s||_{N_\Lambda}$

with

$||(P_\Xi T P_\Xi)^s||_{N_\Lambda}$.

The desired evidence is:

$||T^s||_{N_\Lambda}$ may fail to be contractive,

but

$||(P_\Xi T P_\Xi)^s||_{N_\Lambda}<1$

for $s$ comparable to the effective graph diameter.

This would support the hypothesis that LOD filters out dangerous interface modes.

---

## 10. Optional: Bad-Mode Detection and Small Coarsest Correction

If $|E_0^s|_{G_0}\not<1$, compute dominant generalized singular vectors of $E_0^s$ in the $G_0$ norm.

Solve

$(E_0^s)^*G_0E_0^s x = \lambda G_0x$.

The eigenvectors with $\lambda\ge 1$ identify noncontractive coarse modes.

Store:

* eigenvalue $\lambda$;
* singular value $\sqrt{\lambda}$;
* vector $x$;
* physical LOD function $\Psi x$;
* localization of the mode;
* whether it is boundary-dominated, interface-dominated, or global.

If the number of bad modes is small, construct a small coarsest correction space from these modes and retest $|E_0^s|_{G_0}$ on the complement.

---

## 11. Parameter Sweep

Run the following parameter sweeps.

### Wave numbers

$k\in{16,32,64,128,256}$ if computationally feasible.

### Coarse subdomain diameters

$H_0\in{1,1/2,1/4,1/8}$ or the available equivalents in the code.

### LOD oversampling

$m\in{1,2,3,4,\lceil C\log k\rceil}$.

At minimum, reproduce the experimental setting $m=2$.

### Absorption

$\varepsilon=0$ and $\varepsilon=k$.

### Fine/local scale

Use $\rho\sim k^{-1}$ for the local fine-level decomposition if this is the current experimental setting. (Used in the future, now it is dummy)

---

## 12. Required Output Tables

Produce tables for:

1. Exact two-level GMRES iteration counts.
2. Three-level GMRES iteration counts with one-level coarse preconditioner.
3. Coarse GMRES iteration counts for $A_0$ preconditioned by $M_0^{-1}$.
4. First $s$ such that $|E_0^s|_{G_0}<1$.
5. First $s$ such that $\alpha_s>0$.
6. Values of $|E_0^s|*{G_0}$ for $s=1,\dots,s*{\max}$.
7. Values of $\alpha_s$ for $s=1,\dots,s_{\max}$.
8. Optional: $|(P_\Xi T P_\Xi)^s|*{N*\Lambda}$ versus $|T^s|*{N*\Lambda}$.

---

## 13. Required Plots

Generate plots for each parameter set:

1. GMRES residual history for the exact two-level method.
2. GMRES residual history for the three-level method.
3. $|E_0^s|_{G_0}$ versus $s$.
4. $\alpha_s$ versus $s$.
5. Minimal contraction step $s_0$ versus $H_0^{-1}$.
6. Minimal contraction step $s_0$ versus $k$.
7. Optional: projected and unprojected trace-transfer norms versus $s$.

---

## 14. Main Hypotheses to Verify

The experiments should test the following hypotheses.

### Hypothesis 1: coarse power contraction

There exists $s_0$ such that

$|E_0^{s_0}|_{G_0}<1$.

### Hypothesis 2: geometric iteration scaling

The first contraction step $s_0$ scales mainly like

$s_0\sim \operatorname{diam}(\mathcal G_{\mathrm{eff}})$,

not like $k$.

### Hypothesis 3: polynomial FOV recovery

Even if $S_0=M_0^{-1}A_0$ has no positive field-of-values lower bound, the polynomially transformed operator

$P_s=I-E_0^s$

may satisfy

$\alpha_s>0$

for moderate $s$.

### Hypothesis 4: LOD trace filtering

The projected trace-transfer operator

$P_\Xi T P_\Xi$

is more contractive than the full transfer operator $T$.

### Hypothesis 5: small number of bad modes

If contraction fails, the number of noncontractive singular modes of $E_0^s$ is small enough to motivate a small coarsest correction.

---

## 15. Implementation Warnings

1. Use conjugate transpose for all complex matrices.
2. Use a positive definite norm matrix $G_0$, not the indefinite Helmholtz matrix.
3. Keep left and right preconditioning conventions consistent.
4. Report whether residuals are measured in Euclidean norm, $G_0$ norm, or preconditioned residual norm.
5. Do not treat rank deficiency automatically as an error. It may indicate redundant LOD traces or dependent corrected basis functions. Diagnose and compress.
6. If local matrices are nearly singular, report the condition number and the local domain, wave number, $m$, and $\varepsilon$.
7. In the first stage, prefer exact local solves for $A_{0,\ell}$ so that the behavior of the coarse Schwarz preconditioner is isolated.
8. Save all matrices or function handles needed to reproduce $E_0$, $G_0$, $A_0$, and `apply_M0inv`.

---

## 16. Success Criteria

A successful implementation should answer the following questions.

1. Does the three-level method reproduce the observed $k$-robust GMRES behavior?
2. Does $|E_0^s|_{G_0}$ become less than $1$ for moderate $s$?
3. Is the required $s$ controlled mainly by the number or diameter of coarse subdomains?
4. Does $P_s=I-E_0^s$ recover a positive field-of-values lower bound?
5. Does increasing $m$ improve $|E_0^s|_{G_0}$ or $\alpha_s$?
6. If the power-contraction test fails, are the bad modes low-dimensional?
7. Does adding a small bad-mode coarsest correction improve the contraction?

The most important outputs are:

$s_0=\min{s:|E_0^s|_{G_0}<1}$,

$\alpha_s=\lambda_{\min}\left(\frac{G_0(I-E_0^s)+(I-E_0^s)^*G_0}{2},G_0\right)$,

and GMRES iteration counts as functions of $k$, $H_0$, $m$, and $\varepsilon$.
