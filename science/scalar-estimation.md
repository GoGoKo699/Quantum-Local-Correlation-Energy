# Preparing and estimating the correlation energy

[Gap and signal](gap-and-signal.md) · [Resources](resources-and-regimes.md) · [Home](../README.md)

The output is one energy density. Preparing an entire amplitude vector and then reading every entry would destroy the intended task-level economy.

## An unnormalized source preparation

Start with the maximally correlated cell state

```math
\frac1{\sqrt V}\sum_y|y,y\rangle.
```

Apply the selected $`J/\alpha_J`$ block to the first register. Remove $`x=y`$, sort the two cell labels, and keep an order flag. Since $`J`$ is real symmetric, a Hadamard on that flag coherently adds the two orderings and selects

```math
\frac{b}{\sqrt V B_0},\qquad B_0=\alpha_J/\sqrt2.
```

The remaining norm occupies explicitly flagged failure subspaces. The circuit is unitary and is never conditioned on a measured successful preparation. Thus its cost does not contain an uncharged inverse source-success probability. Small or zero local $`\gamma_x`$ is allowed.

This is a block preparation, not a normalized oracle for an arbitrarily supplied vector. Its cost follows from the specific $`J`$ construction.

## A general positive-gap inverse

Suppose $`A\succeq\Delta I`$ and the full block encoding has normalization $`\alpha_A`$. Standard quantum singular value transformation can approximate

```math
F=\frac{\Delta}{2}A^{-1}
```

on the promised spectral interval using degree

```math
\widetilde O(\alpha_A/\Delta).
```

Here the omitted logarithms include the inverse-approximation accuracy. The source preparation, its inverse, the selected inverse block, and all intervening success projectors give a scalar block

```math
\frac{b^\dagger F b}{VB_0^2}
=-\frac{\Delta e}{2B_0^2}.
```

The corresponding Hadamard-test probability can be estimated coherently by amplitude estimation. For fixed confidence, the leading number of scalar-block invocations is $`O(B_0^2/(\Delta\epsilon))`$. Confidence amplification contributes its own logarithmic factor.

Let $`C_A,C_B`$ be the actual controlled circuit costs of the operator and source constructions. A complete leading-cost statement is

```math
C_{\mathrm{setup}}+
\widetilde O\left[
\frac{B_0^2}{\Delta\epsilon}
\left(C_B+\frac{\alpha_A}{\Delta}C_A\right)
\right].
```

Preparation, inverse preparation, arithmetic uncomputation, scalar estimation and output are included. The formula exposes the two appearances of the inverse gap. It is not a gate count until $`C_A,C_B`$, error allocation and compilation are specified.

The inverse and amplitude-estimation primitives are inherited tools; see [related work](../references/related-work.md).

## A direct finite Neumann implementation

For the constant-on-site example, rank preconditioning supplies

```math
A=D^{1/2}(I-M)D^{1/2},\qquad \|M\|\le q<1.
```

Let $`c=D^{-1/2}b`$. The source is entirely in doubles, so

```math
c=\frac{b}{\sqrt{2E_0}},\qquad
B_c^2=\frac{B_0^2}{2E_0}.
```

Consequently

```math
-e=\frac1V c^\dagger(I-M)^{-1}c.
```

For a degree $`p`$ truncation,

```math
S_p=\sum_{k=0}^p M^k,\qquad
\left|-e-\frac1V c^\dagger S_pc\right|
\le B_c^2\frac{q^{p+1}}{1-q}.
```

This is the standard Neumann-series mechanism. Its usefulness here is that both its physical remainder and its circuit normalization can be bounded from the explicit SD decomposition.

Using the known Green-column envelope $`Z`$, the normalized operator construction has

```math
\beta=
\frac{6(t_o+t_v)}{E_0}
+\frac{3gZ}{E_0}
+\frac{2\alpha_J}{E_0}
+\frac{\alpha_K}{\sqrt2E_0}.
```

Thus $`U_M`$ encodes $`M/\beta`$. Define

```math
\mathcal L_p=\sum_{k=0}^p\beta^k.
```

Preparing a finite LCU label with weights $`\beta^k/\mathcal L_p`$, applying $`k`$ selected $`U_M`$ blocks, and unpreparing the label encodes $`S_p/\mathcal L_p`$. Every product uses fresh success ancillas, so failure components cannot re-enter a good block.

A sequential controlled schedule uses at most $`p`$ operator calls in a SELECT invocation. The coefficient preparation, controls, intermediate projectors, reflections, adjoints and labels are paid in the resource ledger.

The resulting scalar block is

```math
\frac{c^\dagger S_pc}{VB_c^2\mathcal L_p}.
```

Its estimation costs

```math
\widetilde O\left(
\frac{B_c^2\mathcal L_p}{\epsilon}
\,[C_B+pC_M+C_{\mathrm{LCU}}]
\right),
```

plus setup. When $`\beta<1`$, $`\mathcal L_p\le1/(1-\beta)`$. If $`q<1`$ but $`\beta\ge1`$, the physical series still converges while this LCU normalization can become expensive. The example satisfies both inequalities; neither is silently inferred from the other.

## Finite arithmetic and Hermiticity

The coordinate QFT, local input programs, phases, coefficient comparisons, finite LCU states, valid-occupation masks and fermionic signs all have finite precision. A construction that gives a slightly non-Hermitian approximation to a real charged potential can be replaced by the block encoding of $`(\widetilde A+\widetilde A^\dagger)/2`$ at unchanged normalization and constant circuit overhead. Gram products formed with the required selected blocks are Hermitian automatically.

The factor-level errors are bounded in operator norm. Controlled column maps retain their coordinate labels, diagonal multipliers have a uniform error bound, and the source is obtained by contractions of a normalized state. None of these operations requires multiplying the coefficient error by the number of SD matrix entries.

For the Neumann implementation, one may propagate a normalized block error through each finite product, sum with the LCU weights, and add source, coefficient-state, phase-estimation and finite-remainder errors. The [resource ledger](resources-and-regimes.md) gives an explicit conservative allocation rather than treating arithmetic as an exact black box.

## What the scalar output does and does not provide

The algorithm returns an estimate and a stated confidence level for the finite-model LCCSD energy density. It does not return the full excitation vector, validate LCCSD against the exact electronic ground state, or determine the coefficients from a molecular geometry. Those are separate computational tasks and receive separate costs if requested.

The same scalar-output restriction also helps classical algorithms substantially. The [matched comparison](classical-comparison.md) uses it directly.
