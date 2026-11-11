# An explicit LCCSD compiler and its classical boundary

[Home](../README.md) · [Model](model-and-task.md) · [Related work](../references/related-work.md)

## The task and bottleneck

The starting algorithm is the quantum linear-algebra treatment of linearized coupled cluster in Chen and Chan's electronic-structure framework. Its useful output here is a correlation energy per cell, rather than a full classical list of excitation amplitudes.

A fast abstract inverse is insufficient. An implementation must obtain its coefficients, prepare the correlation source, encode the complete SD matrix with a useful normalization, control the inverse gap, and estimate the scalar. A favorable query count can disappear at any of those interfaces.

We address those interfaces for an explicit family of effective electronic models. The primitive interaction is the periodic inverse lattice Laplacian. Each electronic cell contains two primitive sites and one occupied/virtual orbital pair. Their mixing coefficients may vary with position, and the supplied occupied and virtual Fock operators include bounded local hopping.

The input is a finite model with paid coefficient programs. It is not an arbitrary chemistry dataset hidden inside a unit-cost oracle.

## Proposition

Let $L=2^n$, $V=L^3/2$, and let the orbitals, electrostatics and stationary reference satisfy [the model definition](model-and-task.md). Assume bounded local Fock degree and known coefficient bounds, with reversible input-program costs $C_p,C_f$.

Then the following constructions hold.

1. The SD commutator matrix is Hermitian:

   $$
   A=Q_{SD}(H-E_{\mathrm{HF}})Q_{SD}.
   $$

2. Neutral exchange admits normalization

   $$
   \alpha_J=4\pi g\gamma_*^2.
   $$

3. Charged pair potentials admit normalization $gZ$, where the explicit dyadic construction has $Z=7\pi(1-1/L)$. A padded circuit variant has the slightly larger constant $Z=8\pi$.

4. The single-to-double block admits normalization

   $$
   \alpha_K=4g\gamma_*\sqrt{2\pi Z}.
   $$

   Its construction combines a bounded source-column norm with at most two predecessor singles per double.

5. The complete matrix and source therefore admit

   $$
   \alpha_A\le\alpha_F+6gZ+4\alpha_J+\alpha_K,
   \qquad
   B_0=\alpha_J/\sqrt2,\qquad
   \|b\|/\sqrt V\le B_0.
   $$

6. With the certified positive gap $\Delta$ from [the gap proof](gap-and-signal.md), the scalar algorithm has leading cost

   $$
   C_{\mathrm{setup}}+
   \widetilde O\left[
   \frac{B_0^2}{\Delta\epsilon}
   \left(C_B+\frac{\alpha_A}{\Delta}C_A\right)
   \right].
   $$

   Here $C_A,C_B$ include input evaluation, Fourier transforms, local arithmetic, signs, masks, coherent success projections, controls and inverses. Confidence and numerical precision contribute logarithms. No amplitude-vector output cost is hidden.

At fixed physical parameters, fixed positive gap and fixed requested error, every displayed normalization is independent of volume. If $C_p,C_f$ are polynomial in $\log V$ and the bit precision, the complete coherent program has polylogarithmic volume dependence. Setup retains all costs of obtaining and validating those programs.

The proof is distributed across [Coulomb access](coulomb-access.md), [the SD compiler](singles-doubles-compiler.md), [gap and signal](gap-and-signal.md), and [scalar estimation](scalar-estimation.md). The [resource account](resources-and-regimes.md) specializes the inverse to a finite rank-preconditioned Neumann series and charges its arithmetic.

## Why the construction changes the cost

The primitive Coulomb convolution matrix has a norm proportional to $L^2$. Encoding that matrix before using transition-charge neutrality can inherit a size-dependent normalization. The neutral synthesis map instead factors as a discrete divergence, so the two gradients around the inverse Laplacian become an orthogonal projector.

The charged direct interaction is a different operator: a diagonal potential on two retained position labels. Its bounded Green-column overlap does not claim a constant-normalized encoding of the full Coulomb convolution.

The remaining difficulty is the block that adds a local pair to a charge-separated single. Its column may have many nonzero targets, but the column has a bounded squared norm and each target has at most two predecessors. Row/column overlap uses that incidence structure without a factor from the number of cells.

Finally, the source uses the same neutral Gram construction, and the output contracts the inverse with that source coherently. The construction carries the normalization benefit all the way to the requested scalar.

The neutral fixed-rank principle and generic block-encoding tools have prior antecedents. The explicit procedural Green access, SD interface and complete scalar cost are the specific construction assembled and checked here; see [the attribution boundary](../references/related-work.md).

## A finite inverse implementation

For constant on-site separation $E_0$, set $D=mE_0$ at excitation rank $m$. The preconditioned correction

$$
M=I-D^{-1/2}AD^{-1/2}
$$

has a certified physical norm $q<1$ in the example. A separately computed circuit normalization $\beta<1$ makes the finite Neumann polynomial inexpensive to combine by LCU. The physical $q$ controls the truncation error; $\beta$ controls the coefficient sum. This distinction prevents the physical norm from being silently substituted for the actual encoding cost.

The circuit-oriented construction uses uniform-register comparisons for bounded amplitudes and a padded dyadic frequency envelope. The resulting resource figures are conservative analytical circuit budgets, with explicit missing hardware assumptions. They are not measurements of a compiled fault-tolerant implementation.

## Strong classical alternatives change the conclusion

The same weak-coupling gap that makes the quantum inverse well behaved helps a classical variational estimate. Write

$$
F=-e=\frac1Vb^\dagger A^{-1}b.
$$

For any classical trial vector $t$, put $r=b-At$ and

$$
\mathcal F(t)=\frac{2\operatorname{Re}b^\dagger t-t^\dagger At}{V}.
$$

Completing the square gives

$$
0\le F-\mathcal F(t)
=\frac{r^\dagger A^{-1}r}{V}
\le\frac{\|r\|^2}{V\Delta}.
$$

For the constant-on-site example, a diagonal-MP2 trial occupies only determinants with two local pairs. Its full residual nevertheless includes local doubles, singles from $K^\dagger$, and charge-transfer doubles from hopping. All these contributions can be computed from pair matrices and FFTs in $O(V^2\log V)$ work, without storing the unrestricted $O(V^4)$ vector.

The [implemented comparison](classical-comparison.md) reaches the requested nonzero energy tolerance for $V=256$ and $V=2048$. The exact residual identity supplies the mathematical certificate; ordinary double precision supplies the reported evaluations, without a formal interval-rounding guarantee.

Even the diagonal-MP2 scalar alone has an exact $O(V\log V)$ FFT route for arbitrary nonuniform mixing coefficients. Comparing only with a full-vector solver or uncorrected Monte Carlo would therefore overstate the quantum opportunity.

## Scope of the result

The compiler removes a power-of-volume normalization obstacle for this explicit input family and counts the entire scalar algorithm. It does not establish a practical crossover against the strongest classical scalar methods.

The test pattern is a short bit-reversal coefficient program. It is spatially nonuniform but also has inexpensive approximations after discarding low-significance output bits. It is a reproducible algebra and cost fixture, not evidence of hard disorder.

A stronger advantage claim would require a matched regime that survives classical locality, FFT, perturbative, low-rank and stochastic methods; affordable quantum arithmetic and fault tolerance; obtainable model inputs; and a useful total error budget. A fixed-rank norm theorem or a black-box mean-estimation lower bound does not supply those missing facts.

The result applies to SD. The commutator identity used here does not extend automatically to triples, and the finite effective model does not supply a first-principles material-accuracy certificate.
