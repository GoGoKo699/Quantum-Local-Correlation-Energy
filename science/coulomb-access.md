# Coulomb access with explicit normalization

[Model](model-and-task.md) · [SD compiler](singles-doubles-compiler.md) · [Home](../README.md)

The primitive Coulomb matrix has norm growing as $`L^2`$. The construction avoids encoding that whole matrix at constant normalization. It uses two different objects: a bounded neutral response and a controlled diagonal pair potential.

## A uniform diagonal bound

Define

```math
G_0=\frac{4\pi}{N}\sum_{k\ne0}\frac1{\lambda(k)}.
```

Then $`C_{uu}=gG_0`$ and positivity implies $`|C_{uw}|\le gG_0`$. A convenient explicit bound is

```math
G_0\le\overline G=3\pi+\frac{\pi^3}{24}<10.717.
```

To prove it, use centered integer frequencies $`m_j\in[-L/2,L/2)`$ and

```math
\lambda(k)\ge16|m|_2^2/L^2.
```

An $`\ell_\infty`$ shell of radius $`r`$ contains at most $`24r^2+2`$ points. Thus

```math
G_0\le\frac\pi{4L}
\left(24\lfloor L/2\rfloor+2\sum_{r\ge1}r^{-2}\right)
\le3\pi+\frac{\pi^3}{12L}
\le\overline G.
```

The bound is deliberately conservative. Its role is a uniform certificate, not an estimate fitted to finite sizes.

## Neutral exchange as a projector compression

Let $`S`$ embed cell amplitudes into their disjoint internal $`z`$-edges and $`M_\gamma=\operatorname{diag}(\gamma_x)`$. Orient the gradient so the transition-charge synthesis map is $`R=D^\dagger S M_\gamma`$. The neutral Gram matrix is

```math
J=R^\dagger C R
=4\pi g M_\gamma^\dagger S^\dagger\Pi S M_\gamma,
```

where

```math
\Pi=DL_d^+D^\dagger
```

is an orthogonal projector. In momentum space, its nonzero blocks are $`d(k)d(k)^\dagger/\lambda(k)`$, with $`d_j(k)=e^{ik_j}-1`$. Consequently

```math
\|J\|\le\alpha_J=4\pi g\gamma_*^2.
```

The position-dependent $`M_\gamma`$ is a local arithmetic operation. The Fourier transform diagonalizes the primitive electrostatics; it does not require the nonuniform orbital Hamiltonian to be translation invariant.

An equivalent rectangular factor is $`W=C^{1/2}R`$. With the cell-to-site embedding written as $`S_0`$,

```math
W=\sqrt{4\pi g}\,
\operatorname{diag}\!\left(
\frac{e^{-ik_z}-1}{\sqrt{\lambda(k)}}\right)
\mathcal F S_0 M_\gamma.
```

The multiplier has magnitude at most one, so a block encoding uses

```math
\alpha_W=\sqrt{4\pi g}\gamma_*.
```

Local embedding, a coordinate QFT, coefficient arithmetic, and a bounded multiplier implement it. The $`k=0`$ branch is explicitly zero. Its smallest nonzero denominator is inverse polynomial in $`L`$, requiring additional $`O(\log L)`$ working bits at fixed final accuracy.

Constant neutral-Coulomb normalization within a fixed excitation sector is already established in [Chen–Sun–Chan's BSE work](../references/related-work.md). The construction here fixes a procedural input and supplies the particular access needed by the complete SD scalar task.

## Green-half columns from a known envelope

Charged pair potentials need the overlaps of Green-half columns. Directly normalizing such a column by $`G_0`$ would introduce an uncomputed global sum. Instead use a known dyadic envelope.

For signed integer frequencies, let $`B_s=[-2^s,2^s)^3`$, $`S_0=B_0\setminus\{0\}`$, and $`S_s=B_s\setminus B_{s-1}`$ for $`s\ge1`$. Then

```math
|S_s|=7\,8^s,
\qquad s=0,\ldots,n-1,
```

and on shell $`s`$,

```math
\frac{4\pi}{N\lambda(k)}\le\frac\pi{L4^s}.
```

The known squared envelope norm is

```math
Z=\sum_s7\,8^s\frac\pi{L4^s}
=7\pi(1-1/L).
```

Prepare shell $`s`$ with probability $`2^s/(L-1)`$ and a uniform point in that shell. For $`s\ge1`$, 56 subcubes of side $`2^{s-1}`$ enumerate it; $`s=0`$ has seven points. This uses a logarithmic address and finite arithmetic, not an $`N`$-entry frequency table.

A bounded coefficient correction

```math
r_s(k)=\frac{2^{s+1}}{L\sqrt{\lambda(k)}}\le1
```

produces the unnormalized success column

```math
|q_u\rangle=
\sqrt{\frac{4\pi}{NZ}}
\sum_{k\ne0}\frac{e^{-ik\cdot u}}{\sqrt{\lambda(k)}}|k\rangle.
```

Its norm is at most one. The unused norm remains in a failure subspace. No measurement, postselection, amplitude amplification, or knowledge of $`G_0`$ is required.

These columns obey

```math
gZ\langle q_u|q_w\rangle=C_{uw}.
```

A charged orbital density is a coherent combination of its two primitive columns, with weights $`p_x`$ and $`1-p_x`$. Their absolute weights sum to one. Route labels are coherently combined and uncomputed; retaining them as distinguishable garbage would implement a mixture instead of the desired sum.

## The operator encoded by the charged overlap

The particle positions $`u,w`$ remain as control registers. The overlap acts on them as

```math
\sum_{u,w}\frac{C_{uw}}{gZ}
|u,w\rangle\langle u,w|.
```

This is a diagonal pair-potential operator with norm bounded by $`G_0/Z`$. It is not the matrix $`C/(gZ)`$ acting from $`|w\rangle`$ to $`|u\rangle`$. Confusing these operators would incorrectly remove the primitive Coulomb matrix's $`L^2`$ norm.

The density-weighted occupied/virtual pair potentials are obtained by local controlled maps before and after this diagonal operator. Their normalizations remain at most $`gZ`$.

## Success flags are part of the construction

For a block encoding $`P_LUP_R`$, products require intermediate success projections. In particular, the Gram construction is the product of the selected blocks representing $`W^\dagger`$ and $`W`$. The full product $`U_W^\dagger U_W`$ is the identity and is not an encoding of $`W^\dagger W`$.

Likewise, the overlap of two complete Green-column preparation unitaries includes their failure branches unless these branches are removed. Each required success projection is implemented coherently with fresh signal ancillas. All inverse operations, signals and reflections are included in the [scalar resource account](scalar-estimation.md).

The reciprocal coefficient functions, phase factors, masks and uncomputation have explicit arithmetic cost. The [resource page](resources-and-regimes.md) gives a circuit-oriented variant and quantitative budgets for these operations.
