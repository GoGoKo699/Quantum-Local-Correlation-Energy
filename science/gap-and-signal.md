# A certified inverse gap and a nonzero signal

[SD compiler](singles-doubles-compiler.md) · [Scalar estimation](scalar-estimation.md) · [Home](../README.md)

An inexpensive block encoding is insufficient if the inverse is nearly singular or if the requested accuracy permits returning zero. This page certifies both issues for an explicit parameter regime.

## A general sufficient gap

Let

```math
E_g=\lambda_{\min}(f_v)-\lambda_{\max}(f_o)>0.
```

The Fock excitation operator is at least $`mE_g`$ on rank $`m`$, hence at least $`E_g`$ throughout SD. The [interaction decomposition](singles-doubles-compiler.md) gives

```math
\|W_{SD}\|\le\Gamma,\qquad
\Gamma=g\left[
6\overline G+16\pi\gamma_*^2
+4\gamma_*\sqrt{2\pi\overline G}\right].
```

Therefore

```math
A\succeq\Delta I,\qquad
\Delta=E_g-\Gamma>0.
```

This is a sufficient condition, not a claim that every gapped mean-field model satisfies it. A computable Fock lower bound follows from diagonal bounds and hopping row sums:

```math
E_g\ge
\min_x(f_v)_{xx}-\max_x(f_o)_{xx}
-t_{v,\mathrm{row}}-t_{o,\mathrm{row}}.
```

The global bounds are supplied with, or proved from, the input program. Discovering them by searching arbitrary input data is a separate setup cost.

## A stronger rank-aware bound for constant on-site separation

Suppose the occupied and virtual on-site energies are $`-E_0/2`$ and $`+E_0/2`$, with at most six neighbors per band and hopping magnitudes at most $`t_o,t_v`$. Define $`D=mE_0`$ on excitation rank $`m=1,2`$ and write

```math
A=D+T+W,\qquad
M=-D^{-1/2}(T+W)D^{-1/2}.
```

Let

```math
\overline K=4g\gamma_*\sqrt{2\pi\overline G}.
```

Dividing each block by its actual rank denominator yields

```math
\|M\|\le q=
\frac{6(t_o+t_v)}{E_0}
+\frac{3g\overline G}{E_0}
+\frac{2\alpha_J}{E_0}
+\frac{\overline K}{\sqrt2E_0}.
```

The direct term has at most six pairs in doubles, divided by $`2E_0`$; the single-sector bound is smaller. Four exchange lifts in doubles likewise give $`2\alpha_J/E_0`$. The off-diagonal $`K`$ block is divided by $`\sqrt{E_0(2E_0)}`$. The hopping contribution has $`m`$ particles of each type and its factor $`m`$ cancels the denominator.

If $`q<1`$, then

```math
A=D^{1/2}(I-M)D^{1/2}
\succeq(1-q)D
\succeq\Delta_{\mathrm{rank}} I,
\qquad
\Delta_{\mathrm{rank}}=E_0(1-q).
```

This bound benefits both the quantum inverse and the classical residual certificate. Replacing $`\overline G`$ by a paid, certified finite-volume $`G_0`$ can sharpen it. Such a replacement is optional; the uniform bound already suffices for the numerical regime below.

The physical contraction $`q`$ differs from the circuit normalization $`\beta`$ in [scalar estimation](scalar-estimation.md). A small physical norm does not automatically make an arbitrary block encoding inexpensive.

## Source norm without a volume-dependent preparation probability

Only doubles containing two local electron–hole pairs have nonzero source coefficients. With consistent determinant phases,

```math
b_{xy}=J_{xy},\qquad x<y.
```

The bounded Gram encoding supplies the source envelope

```math
\frac{\|b\|^2}{V}\le B_0^2,\qquad
B_0=\frac{\alpha_J}{\sqrt2}.
```

No lower bound on $`\|b\|`$ is needed by the preparation or scalar estimator. A lower bound is useful only to certify that the example's output task is nontrivial.

## A uniform lower signal for a restricted example family

In this subsection, take $`L\ge8`$ and $`1/4\le p_x\le3/4`$. Then $`\gamma_{\min}=\sqrt3/4`$ and $`\gamma_*=1/2`$. Write

```math
J=4\pi g M_\gamma K_0 M_\gamma.
```

The primitive electrostatics makes $`K_0`$ circulant on the cell grid, even though nonuniform $`M_\gamma`$ makes $`J`$ noncirculant. Its eigenvalue, obtained by aliasing the two primitive $`z`$ frequencies, is

```math
j(k_x,k_y,\theta)=\frac12
\left[\frac{u}{a+u}+\frac{1-u}{a+1-u}\right],
```

```math
a=\sin^2(k_x/2)+\sin^2(k_y/2),\qquad
u=\sin^2(\theta/4).
```

A zero-denominator fraction is set to zero, as prescribed by the removed primitive zero mode.

At least $`1/32`$ of the finite cell frequencies have $`|k_x|,|k_y|\le\pi/4`$ and $`\theta\in[\pi/2,3\pi/2]`$. On this set $`j\ge8/15`$. At least $`1/4`$ have $`|k_x|,|k_y|\ge\pi/2`$, on which $`j\le1/3`$. These fractions hold for powers of two $`L\ge8`$.

For the first bound, set $`a_0=1-1/\sqrt2`$. The displayed function decreases with $`a`$ and is concave and symmetric in $`u`$, so its minimum on the first set occurs at $`a=a_0`$, $`u=a_0/2`$ or $`1-a_0/2`$, and exceeds $`8/15`$. On the second set, $`a\ge1`$ and the maximum over $`u`$ is $`1/(2a+1)\le1/3`$.

If a random variable occupies two intervals separated by $`d`$ with probabilities at least $`p_1,p_2`$, minimizing its squared deviations over a candidate mean gives variance at least $`p_1p_2d^2/(p_1+p_2)`$. Hence

```math
\mathrm{Var}(j)\ge
\frac{(1/32)(1/4)}{1/32+1/4}\left(\frac15\right)^2
=\frac1{900}.
```

Circulant Parseval and the constant diagonal give

```math
\frac1V\sum_{x\ne y}|(K_0)_{xy}|^2
=\mathrm{Var}(j).
```

Using the unordered source and $`\gamma_{\min}\le\gamma_x`$,

```math
\frac{B_0^2}{1600}\le\frac{\|b\|^2}{V}\le B_0^2.
```

This lower bound holds for every spatial arrangement of the allowed $`p_x`$, including a deterministic nonperiodic pattern.

If $`A\preceq\Lambda I`$, positivity therefore gives

```math
\frac{B_0^2}{1600\Lambda}\le |e|\le\frac{B_0^2}{\Delta}.
```

One conservative upper bound is $`\Lambda=\alpha_F+\Gamma`$. The stronger rank-aware gap may replace $`\Delta`$ on the right.

## Numerical regime used throughout

Use $`g=0.01`$, $`E_0=2`$, $`t_o=t_v=0.02`$, and $`\gamma_*=1/2`$. Then

```math
E_g\ge1.76,\quad
\Delta\ge0.827218,\quad
\Lambda\le5.412782,\quad
B_0^2=4.9348022\times10^{-4}.
```

The rank-aware gap is larger; its unrounded value is produced by the [resource calculator](../verification/README.md). Set

```math
\epsilon=\frac{B_0^2}{16000\Lambda}
=5.6980891\times10^{-9}.
```

For $`L\ge8`$, the true correlation energy magnitude is at least ten times this accuracy. This rules out a trivial zero answer. The accuracy is in the effective model's energy units, not a declared material-level chemical-accuracy target.

The [classical comparison](classical-comparison.md) shows that this nonzero task nevertheless has an inexpensive classical solution. A nonzero signal and a certified gap do not establish quantum advantage.
