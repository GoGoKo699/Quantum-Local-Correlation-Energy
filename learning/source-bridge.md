# From scientific computing to correlation energy

[Home](../README.md) · [Milestone map](milestones.md) · [Scientific overview](../science/overview.md)

## The one educational anchor

Lin Lin, *Lecture Notes on Quantum Algorithms for Scientific Computation*, **arXiv:2201.08309v1**: [version record](https://arxiv.org/abs/2201.08309v1), [complete PDF](https://arxiv.org/pdf/2201.08309v1). The submission is dated 20 January 2022; the PDF title page is dated 21 January 2022. The pointers below use this exact version and its printed page numbers. Later course-note revisions are not substituted.

The computational question is concrete: how much does allowing coupled electronic excitations lower the energy per cell in the supplied finite model? The project expresses that correction as

```math
e=-\frac{b^\dagger A^{-1}b}{V}.
```

Here $`b`$ couples the reference to double excitations, $`A`$ is the singles-and-doubles response matrix, and $`V`$ is the number of electronic cells. The output is one scalar with a requested additive error. [The model definition](../science/model-and-task.md) specifies the inputs, stationary reference, and distinction from exact material energies.

The anchor supplies the quantum linear-algebra vocabulary. The project must still construct the particular matrix and source, control their errors, and compare complete ways of obtaining the same scalar. Its supported result is an explicit compiler for a restricted model family, together with a strong classical countercheck. It establishes no practical quantum crossover.

## A targeted route through the anchor

These are entry points for the [30 project milestones](milestones.md), not a requirement to read the notes in full. Section numbers are the stable locators; the page column gives their starting pages.

| Anchor sections | Pages | Use in this project | Canonical project account |
|---|---:|---|---|
| §§1.1, 1.3, 1.5 | 7, 13, 16 | Interpret amplitudes, registers, selected branches, and measured probabilities. | [Source preparation](../science/scalar-estimation.md#an-unnormalized-source-preparation) |
| §§1.6–1.10 | 19–24 | Account for finite arithmetic, uncomputation, accumulated error, and the difference between queries and a circuit budget. | [Resource account](../science/resources-and-regimes.md) |
| §§3.1, 4.2; Example 4.5 | 39, 57–58 | Convert a controlled overlap into the energy estimate and charge coherent repetitions. | [Scalar estimation](../science/scalar-estimation.md) |
| §§3.3, 4.4 | 45, 65 | Connect Fourier coordinates and inverse differential operators to the project's specified periodic lattice kernel. | [Coulomb access](../science/coulomb-access.md) |
| §§6.1–6.3; Definition 6.1 | 79–84 | Identify what an access routine actually implements, including its scale and approximation error. | [Paid input model](../science/model-and-task.md#paid-access-and-error-contract) |
| §6.4, §7.3; Lemma 7.9 | 85, 99 | Combine sectors and a finite inverse polynomial while retaining success flags. | [SD compiler](../science/singles-doubles-compiler.md), [finite inverse](../science/scalar-estimation.md#a-direct-finite-neumann-implementation) |
| §§7.5–7.6, 8.3–8.4 | 105, 109, 122, 124 | Understand the general polynomial-inverse route and its spectral requirements. | [Positive-gap inverse](../science/scalar-estimation.md#a-general-positive-gap-inverse) |
| §4.3; Remark 4.6 | 59–61 | Separate a normalized solution state from its norm and from the classical number actually requested here. | [Output contract](../science/scalar-estimation.md#what-the-scalar-output-does-and-does-not-provide) |

## Translate the notation before importing a bound

| Anchor convention | Project translation | Consequence |
|---|---|---|
| An $`n`$-qubit matrix acts on a space of dimension $`2^n`$. | $`N=L^3`$ counts primitive sites, $`V=N/2`$ counts cells, and the SD determinant space has another dimension. | A logarithmic number of address qubits does not make acquisition of arbitrary coefficient data inexpensive. |
| A block encoding has a normalization $`\alpha`$ and an operator error. | $`U_A`$ encodes $`A/\alpha_A`$; $`\alpha_J`$, $`gZ`$, and $`\alpha_K`$ come from different constructions. | An upper bound on the physical norm need not be an achievable circuit normalization. |
| The inverse construction uses a spectral cutoff for a scaled matrix. | The relevant cutoff is $`\Delta/\alpha_A`$ when $`A\succeq\Delta I`$. | A one-particle orbital separation is not by itself the required SD inverse gap. |
| The linear-system right-hand side is often a normalized state. | The good branch here is $`b/(\sqrt V B_0)`$, with explicitly retained failure branches. | The algorithm need not know $`\lVert b\rVert`$ or postselect successful source preparation. |
| Section 8.4 uses $`\beta`$ as a reciprocal-polynomial scaling factor. | This project's $`\beta`$ normalizes the preconditioned correction $`M`$. | These are different quantities. Here $`q`$ bounds $`\lVert M\rVert`$, while $`\beta`$ controls the implemented LCU coefficient sum. |
| A linear solver may return a state proportional to $`A^{-1}b`$. | This task returns $`e`$ and a confidence statement. | Neither a full amplitude list nor a material-accuracy validation is included in that output. |

## A one-component view of the scalar mechanism

Consider a single response coordinate with $`A=a\ge\Delta>0`$, real source $`b`$, a known positive envelope $`B_0`$, and $`|b|\le\sqrt V B_0`$. This is an algebraic illustration of the readout, not an electronic benchmark. A unitary source preparation has selected amplitude $`b/(\sqrt V B_0)`$; a selected inverse block has value $`\Delta/(2a)`$. Composing them with the required intermediate projectors gives the real scalar

```math
s=\frac{b^2}{VB_0^2}\frac{\Delta}{2a}
  =-\frac{\Delta e}{2B_0^2}.
```

A Hadamard test of the resulting unitary embedding has probability $`(1+s)/2`$ for its chosen outcome. Estimating that probability recovers $`s`$, then $`e=-2B_0^2s/\Delta`$. The measured information is sufficient for the requested energy; there is no need to extract the response amplitude separately. The construction remains meaningful at $`b=0`$ because it never divides by the source norm.

Every repetition still uses source preparation, inverse preparation, the inverse block, controls, and coherent success selection. Removing a success projector can change the encoded operator: an unqualified $`U^\dagger U`$ is the identity, not the desired Gram product. The full implementation and its cost are in [scalar estimation](../science/scalar-estimation.md).

## The steps supplied by this project

**A physical response matrix.** The [SD projection argument](../science/model-and-task.md#why-sd-is-hermitian-here) establishes Hermiticity from the supplied stationary reference and two-body Hamiltonian. This is a statement about the chosen excitation space. The corresponding commutator argument does not extend automatically to unrestricted triples.

**An explicit long-range interaction.** The [Coulomb construction](../science/coulomb-access.md) uses the full periodic inverse lattice Laplacian with its zero mode removed. Neutral transition charges turn an inverse-Laplacian sandwich into a projector compression. Charged direct terms use a different, diagonal pair-potential operator with both position labels retained. A bounded encoding of that operator is not a bounded encoding of the full Coulomb convolution matrix. The chosen lattice kernel also does not follow merely by citing the anchor's Poisson example.

**The complete SD interface.** The [compiler](../science/singles-doubles-compiler.md) includes the single-to-double block, not only the rank-preserving interaction. Its column norm and bounded target-predecessor count supply the row/column construction. Distinct failure labels prevent unwanted failure–failure overlaps.

**Two separate spectral certificates.** The [gap and signal analysis](../science/gap-and-signal.md) supplies a sufficient positive SD gap and, for a narrower example family, a lower bound on the energy signal. The first makes inversion stable; the second ensures the chosen absolute-error target cannot be met by reporting zero. Neither certificate establishes useful material accuracy.

**A finite implementation with paid inputs.** The [resource account](../science/resources-and-regimes.md) replaces an abstract inverse in the worked regime with rank preconditioning and a Neumann polynomial of degree 11. Both $`q<1`$ and $`\beta<1`$ are checked separately. Reversible coefficient programs, amplitude comparisons, Fourier operations, projections, precision and confidence enter the budget. The reported gate figures are conservative analytical construction budgets; they are not a compiled fault-tolerant netlist or measured runtime.

**A matched classical scalar comparison.** The [classical calculation](../science/classical-comparison.md) uses a local-pair trial and a full residual certificate, including singles and charge-transfer doubles. It meets the stated error tolerance on the published examples without storing the unrestricted SD vector. The physical gap helps this certificate too. Pair FFTs, locality, and the short coefficient program must therefore be considered before interpreting favorable quantum volume scaling as an advantage. The reported floating-point residual evaluation does not provide a formal interval-rounding guarantee.

## Read the contribution at its supported scope

The complete story joins constant normalization, an unnormalized source, a controlled inverse, and scalar readout under explicit input and gap assumptions. At fixed physical parameters, fixed positive gap and fixed requested error, the construction avoids a power-of-volume normalization cost. The remaining comparison is determined by arithmetic, setup, precision, model usefulness, and the best classical treatment of the same output. In the demonstrated regime the classical method is already effective.

The anchor supports the shared mathematical language. [Related work](../references/related-work.md) attributes the existing quantum LCC framework, earlier neutral fixed-rank ideas, algorithmic primitives and classical methods. Those research citations are proof and attribution dependencies, not additional educational anchors. The [milestone map](milestones.md) connects these distinctions into the project's two-phase learning route.
