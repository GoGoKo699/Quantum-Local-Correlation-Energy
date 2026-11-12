# Model and scalar task

[Scientific overview](overview.md) · [Coulomb access](coulomb-access.md) · [Home](../README.md)

## The quantity to compute

Linearized coupled cluster reduces a correlation calculation to a linear equation in excitation amplitudes. We study singles and doubles, with a supplied stationary Hartree–Fock reference $`|\Phi\rangle`$. For excitation operators $`O_\nu`$ and determinants $`|\nu\rangle=O_\nu|\Phi\rangle`$, define

```math
A_{\mu\nu}=\langle\mu|[H,O_\nu]|\Phi\rangle,
\qquad b_\mu=\langle\mu|H|\Phi\rangle,
\qquad At=-b.
```

The requested output is the correlation energy per electronic cell,

```math
e=\frac{b^\dagger t}{V}
  =-\frac{b^\dagger A^{-1}b}{V},
```

to additive error $`\epsilon`$ and a specified failure probability. Only this scalar is returned. Reading the full vector $`t`$ would be a different output problem.

This is the LCCSD task within the [Chen–Chan framework](../references/related-work.md). LCCSD is an approximate electronic-structure model. Accurate evaluation of its scalar is distinct from accuracy relative to the exact many-electron ground-state energy.

## An explicit input family

The primitive spatial grid is the cubic torus

```math
\Lambda=(\mathbb Z/L\mathbb Z)^3,
\qquad L=2^n,
\qquad N=L^3.
```

The lattice spacing is fixed. Pair $`(x_1,x_2,2x_3)`$ with $`(x_1,x_2,2x_3+1)`$ to form one electronic cell. There are $`V=N/2`$ cells, one spinless occupied orbital and one virtual orbital per cell, and $`V`$ reference electrons.

For a parameter $`0\le p_x\le1`$, the orthonormal local orbitals are

```math
|o_x\rangle=\sqrt{p_x}|x,+\rangle+
\sqrt{1-p_x}|x,-\rangle,
```

```math
|v_x\rangle=-\sqrt{1-p_x}|x,+\rangle+
\sqrt{p_x}|x,-\rangle.
```

Different cells have disjoint primitive support. Set $`\gamma_x=\sqrt{p_x(1-p_x)}`$ and supply a bound $`\gamma_x\le\gamma_*\le1/2`$. The densities are

```math
\rho_{o_x}=p_x\delta_{x,+}+(1-p_x)\delta_{x,-},
\qquad
\rho_{v_x}=(1-p_x)\delta_{x,+}+p_x\delta_{x,-},
```

```math
\rho_x=\rho_{v_xo_x}
=\gamma_x(\delta_{x,-}-\delta_{x,+}).
```

The transition density has exactly zero total charge. The same computed $`\gamma_x`$ multiplies both signs, so a finite numerical representation need not introduce a spurious monopole. No localization-tail approximation is involved in this family.

## Full periodic electrostatics

Let $`D`$ be the forward discrete gradient and $`L_d=D^\dagger D`$ the positive discrete Laplacian. The interaction is

```math
C=4\pi gL_d^+,
\qquad g>0,
```

where the pseudoinverse removes the zero Fourier mode. For $`k_j=2\pi m_j/L`$,

```math
\lambda(k)=4\sum_{j=1}^3\sin^2(k_j/2),
\qquad
\widehat C(k)=
\begin{cases}
4\pi g/\lambda(k),&k\ne0,\\
0,&k=0.
\end{cases}
```

This defines the ultraviolet resolution and periodic zero-mode convention as part of the finite input. It is not the minimum-image function $`1/r`$ and is not a silently truncated continuum integral.

The two-body coefficients are Coulomb inner products of orbital-product densities,

```math
(pq|rs)=\sum_{u,w\in\Lambda}
\phi_p(u)^*\phi_q(u)C_{uw}
\phi_r(w)^*\phi_s(w).
```

The normal-ordered Fock operator is separately supplied as Hermitian occupied and virtual blocks $`f_o,f_v`$, with the cross block zero. Each block has bounded local degree. Its diagonal and hopping coefficients may depend on position. Nonzero hopping allows a virtual electron and its hole to separate, so the SD space includes charge-transfer configurations.

Equivalently, one can choose the one-body matrix of the underlying two-body Hamiltonian so that these are its exact reference Fock blocks. This establishes an electronic model and a stationary reference. It does not establish that an arbitrary prescribed Fock matrix arises from a local nuclear potential, or that deriving it from a material is inexpensive.

## Paid access and error contract

A reversible program evaluates $`p_x`$ at cost $`C_p(b,\log V)`$, and another evaluates a specified Fock neighbor and coefficient at cost $`C_f(b,\log V)`$. The requested bit precision is $`b`$. Known global coefficient bounds follow from the input specification; optimizing arbitrary functions to discover such bounds is not granted for free.

The principal compact-input regime has $`C_p,C_f`$ polynomial in $`b`$ and $`\log V`$. The same programs are available to a classical algorithm. Independently tabulated data instead pay their acquisition and coherent lookup costs. No QRAM is assumed.

The setup ledger includes model construction, reference validation, coefficient-program construction, classical polynomial/LCU compilation, and any reusable data processing. The runtime ledger includes forward and inverse preparations, controlled calls, numerical arithmetic, signs, masks, and scalar estimation. If this effective model is used to approximate a material, its physical-model error receives a separate budget.

## Why SD is Hermitian here

For a two-body Hamiltonian and a stationary reference,

```math
H|\Phi\rangle=E_{\mathrm{HF}}|\Phi\rangle+
\sum_{\eta\in D}b_\eta|\eta\rangle.
```

An excitation operator of positive rank acting on a double either gives zero or an excitation of rank at least three. Therefore the second term is killed by projection onto singles and doubles. With $`Q_{SD}`$ the SD projector,

```math
A=Q_{SD}(H-E_{\mathrm{HF}})Q_{SD}=A^\dagger.
```

This identity supplies the Hermitian scalar problem used throughout the repository. Positivity requires the separate [gap certificate](gap-and-signal.md); a bare mean-field band gap alone is insufficient.

The argument fails for a general larger excitation truncation. A triples-to-singles de-excitation can survive while its reverse disconnected creation is removed by the commutator. Thus neither Hermiticity nor a volume-independent inverse gap is transferred to unrestricted LCCSDT.

## Finite-precision perturbations

Let $`\|b\|/\sqrt V\le B_0`$, $`A\succeq\Delta I`$, and let $`\widetilde A=A+\delta A`$, $`\widetilde b=b+\delta b`$, with

```math
\|\delta A\|\le\eta_A\le\Delta/2,
\qquad \|\delta b\|/\sqrt V\le\eta_b.
```

The resolvent identity and $`\|\widetilde A^{-1}\|\le2/\Delta`$ give

```math
|\widetilde e-e|
\le\frac{2B_0^2\eta_A}{\Delta^2}
 +\frac{4B_0\eta_b+2\eta_b^2}{\Delta}.
```

This is an operator/state-norm error condition. The explicit Fourier construction controls these norms by composition, rather than accumulating a worst-case error over all $`V^4`$ SD entries. Additional error from the inverse approximation, coherent composition and scalar estimator is allocated separately.

The next page constructs the required [Coulomb access](coulomb-access.md).
