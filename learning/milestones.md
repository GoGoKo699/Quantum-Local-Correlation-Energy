# Thirty milestones from scientific computing to the complete scalar task

[Home](../README.md) · [Source bridge](source-bridge.md) · [Scientific overview](../science/overview.md) · [Verification](../verification/README.md)

This map explains how an explicit quantum compiler for a finite electronic model becomes a complete correlation-energy calculation, and why the worked regime still has an inexpensive classical solution. It has exactly 30 understanding targets. They describe the project for any reader; they are independent of a reader's previous study or progress.

The single educational anchor is **Lin Lin, [*Lecture Notes on Quantum Algorithms for Scientific Computation*, arXiv:2201.08309v1](https://arxiv.org/abs/2201.08309v1)**, submitted 20 January 2022. All “Lin” section numbers below refer to that version. The [source bridge](source-bridge.md) supplies the section guide, notation translation and missing steps. Lin provides the computational foundations; the linked project pages supply the electronic model, special compiler proofs and matched classical analysis. Research attribution remains in [related work](../references/related-work.md).

The allocation is **12 background/preliminary milestones and 18 narrative/innovation milestones**. The first phase joins local physical objects, matrix operations and access costs into a fair scalar-computation problem. The second assembles the model-specific compiler, then its error and resource accounts, and finally the classical countercheck and contribution boundary. This allocation gives the interfaces and comparison as much attention as the inverse itself. It does not attempt to cover every lemma or all of quantum chemistry.

## Map at a glance

| Background and preliminaries | Narrative and innovation |
|---|---|
| [M01. Scalar output and accuracy](#m01-scalar-output-and-accuracy) | [M13. The exact SD reduction](#m13-the-exact-sd-reduction) |
| [M02. Local orbitals and neutral transitions](#m02-local-orbitals-and-neutral-transitions) | [M14. Neutrality as a projector compression](#m14-neutrality-as-a-projector-compression) |
| [M03. Excitation labels and fermionic signs](#m03-excitation-labels-and-fermionic-signs) | [M15. Green columns from a known envelope](#m15-green-columns-from-a-known-envelope) |
| [M04. Hermiticity, gaps and inverse sensitivity](#m04-hermiticity-gaps-and-inverse-sensitivity) | [M16. Charged overlaps and their operator](#m16-charged-overlaps-and-their-operator) |
| [M05. Obtainable inputs and reversible arithmetic](#m05-obtainable-inputs-and-reversible-arithmetic) | [M17. The single-to-double incidence bound](#m17-the-single-to-double-incidence-bound) |
| [M06. Periodic kernels and Fourier coordinates](#m06-periodic-kernels-and-fourier-coordinates) | [M18. A coherent compiler for all SD sectors](#m18-a-coherent-compiler-for-all-sd-sectors) |
| [M07. What a block encoding promises](#m07-what-a-block-encoding-promises) | [M19. A source that stays unnormalized](#m19-a-source-that-stays-unnormalized) |
| [M08. Coherent products and success registers](#m08-coherent-products-and-success-registers) | [M20. A gap certified from the supplied model](#m20-a-gap-certified-from-the-supplied-model) |
| [M09. From an overlap to a measured scalar](#m09-from-an-overlap-to-a-measured-scalar) | [M21. Accuracy that excludes the zero answer](#m21-accuracy-that-excludes-the-zero-answer) |
| [M10. Polynomial inverses and normalization](#m10-polynomial-inverses-and-normalization) | [M22. Two bounds for one Neumann inverse](#m22-two-bounds-for-one-neumann-inverse) |
| [M11. Classical scalar certificates](#m11-classical-scalar-certificates) | [M23. Finite amplitudes from comparisons](#m23-finite-amplitudes-from-comparisons) |
| [M12. A complete and matched cost question](#m12-a-complete-and-matched-cost-question) | [M24. One coherent error budget](#m24-one-coherent-error-budget) |
| | [M25. Paid output precision and confidence](#m25-paid-output-precision-and-confidence) |
| | [M26. What the concrete resource ledger establishes](#m26-what-the-concrete-resource-ledger-establishes) |
| | [M27. The classical FFT shortcut](#m27-the-classical-fft-shortcut) |
| | [M28. The full residual in a pair representation](#m28-the-full-residual-in-a-pair-representation) |
| | [M29. Reading the matched numerical evidence](#m29-reading-the-matched-numerical-evidence) |
| | [M30. The contribution and the remaining advantage gap](#m30-the-contribution-and-the-remaining-advantage-gap) |

## Phase I: background and preliminaries

### M01. Scalar output and accuracy

**Phase:** Background/preliminaries.

**Central question:** What useful number is the computation required to return?

**Understanding target:** Explain why the task is the LCCSD correlation energy per cell, with an additive error and failure probability. Distinguish evaluating this finite-model scalar from recovering all excitation amplitudes or predicting a material's exact ground-state energy.

**Dependencies and story role:** No earlier milestone is required. This fixes the output contract before introducing either computational method.

**Conceptual connection:** A scientific quantity becomes an algorithmic task only after its input, accuracy and observable output are specified. A quantum solution state alone is not the requested classical number.

**Sources:** Lin §§1.5, 1.10 and 4.3; [Model: the quantity to compute](../science/model-and-task.md#the-quantity-to-compute); [Scalar output scope](../science/scalar-estimation.md#what-the-scalar-output-does-and-does-not-provide).

### M02. Local orbitals and neutral transitions

**Phase:** Background/preliminaries.

**Central question:** How can two charged orbital densities produce a neutral transition density?

**Understanding target:** Read the occupied and virtual orbitals on a two-site cell, their orthogonality and their density weights. Explain why the transition has equal and opposite weights, even when the mixing parameter varies across cells.

**Dependencies and story role:** Builds on M01 by identifying the local physical objects from which the input is made.

**Conceptual connection:** Inner products and superpositions from quantum mechanics meet a concrete orbital representation. Exact local neutrality later removes a normalization penalty without assuming a translation-invariant electronic Hamiltonian.

**Sources:** Lin §§1.1 and 1.11 for state/vector notation; [Model: explicit input family](../science/model-and-task.md#an-explicit-input-family). The orbital construction is a project supplement.

### M03. Excitation labels and fermionic signs

**Phase:** Background/preliminaries.

**Central question:** What does a single or double excitation record?

**Understanding target:** Describe electron and hole labels relative to a stationary reference, distinguish a local pair from a charge-separated excitation, and explain why sorting and collision masks must preserve determinant signs. Understand that logarithmic label registers do not make the full classical amplitude list small.

**Dependencies and story role:** Uses M02 to move from one cell to the restricted many-electron space.

**Conceptual connection:** Register encodings and reversible label operations implement a physically constrained basis. Fixed excitation rank bounds the number of coordinate registers and routes, while the number of basis states still grows rapidly.

**Sources:** Lin §§1.3 and 1.7 for circuit/register foundations; [Model: SD target](../science/model-and-task.md#the-quantity-to-compute); [Compiler: normal ordering and quasiparticle count](../science/singles-doubles-compiler.md#normal-ordering-and-the-finite-quasiparticle-count). Fermionic details are supplied locally.

### M04. Hermiticity, gaps and inverse sensitivity

**Phase:** Background/preliminaries.

**Central question:** Which matrix properties make a stable inverse-energy calculation possible?

**Understanding target:** Distinguish Hermiticity, positive definiteness, a lower spectral gap and an encoding normalization. Explain why perturbations become more consequential near a small inverse gap, and why a mean-field band gap does not by itself certify the interacting SD matrix.

**Dependencies and story role:** Uses M01 and M03 to identify the linear-algebra properties the physical reduction must prove.

**Conceptual connection:** The spectral viewpoint behind quantum linear solvers becomes a quantitative error contract for the actual electronic matrix.

**Sources:** Lin §§4.3, 6.4 and 8.4; [Model: Hermiticity](../science/model-and-task.md#why-sd-is-hermitian-here); [Model: finite-precision perturbations](../science/model-and-task.md#finite-precision-perturbations); [Gap certificate](../science/gap-and-signal.md#a-general-sufficient-gap).

### M05. Obtainable inputs and reversible arithmetic

**Phase:** Background/preliminaries.

**Central question:** What does it cost to obtain a coefficient in a coherent computation?

**Understanding target:** Distinguish a short reversible coefficient program from a large supplied table. Explain evaluation, retained outputs, uncomputation, precision and controlled calls, and identify model/reference construction as setup work. Give the classical comparator the same input program.

**Dependencies and story role:** Uses M02–M04 to connect physical coefficients to a paid computational interface.

**Conceptual connection:** Lin's query model is an abstraction whose implementation must be charged. A compact program can support favorable scaling; an arbitrary dataset does not inherit that access cost.

**Sources:** Lin §§1.7, 1.8, 1.10 and 6.1; [Paid access and error contract](../science/model-and-task.md#paid-access-and-error-contract); [Resource task and input](../science/resources-and-regimes.md#1-task-input-and-what-is-charged).

### M06. Periodic kernels and Fourier coordinates

**Phase:** Background/preliminaries.

**Central question:** Which part of this nonuniform electronic problem remains a convolution?

**Understanding target:** Explain the finite periodic inverse lattice Laplacian, the removed zero mode, and the distinction between the primitive electrostatic kernel and the nonuniform orbital factors. Recognize that both QFT and classical FFT methods can use the kernel structure.

**Dependencies and story role:** Uses M02 and M05 to turn the spatial input into structured linear algebra.

**Conceptual connection:** Lin's Poisson example motivates spectral access, while the project fixes its own lattice, zero-mode convention and scalar task. No continuum or distance-cutoff replacement is implicit.

**Sources:** Lin §§3.3 and 4.4; [Full periodic electrostatics](../science/model-and-task.md#full-periodic-electrostatics); [Classical local-pair and FFT representation](../science/classical-comparison.md#3-local-pair-trial-and-an-exact-mp2-fft).

### M07. What a block encoding promises

**Phase:** Background/preliminaries.

**Central question:** Why does the normalization of an encoded matrix affect the whole algorithm?

**Understanding target:** Identify the matrix selected by signal projectors from a unitary, its normalization and its approximation error. Distinguish an operator norm bound from the normalization achieved by a particular circuit. Explain how an artificial volume factor can propagate into inverse and scalar-estimation costs.

**Dependencies and story role:** Uses M04–M06 to define the interface the compiler must implement.

**Conceptual connection:** A block encoding is more informative than an entry oracle: its selected operator, normalization, cost and error must all be specified together.

**Sources:** Lin §§6.1–6.4, especially Definition 6.1; [Overview proposition](../science/overview.md#proposition); [Coulomb access](../science/coulomb-access.md).

### M08. Coherent products and success registers

**Phase:** Background/preliminaries.

**Central question:** Why is multiplying preparation unitaries different from multiplying their selected blocks?

**Understanding target:** Explain the intermediate projection needed in a Gram product, why an unconstrained unitary followed by its inverse gives the identity, and why distinct failure sectors or fresh signal registers matter. Distinguish coherent selection from measured postselection.

**Dependencies and story role:** Builds on M07 and M05 to make composition of the forthcoming access primitives trustworthy.

**Conceptual connection:** The inner-product and LCU constructions in the anchor require careful bookkeeping of which subspace represents the intended matrix. Uncomputation removes workspace only when it does not erase necessary success information.

**Sources:** Lin §6.2, including Eq. (6.16), and §7.3; [Success flags](../science/coulomb-access.md#success-flags-are-part-of-the-construction); [Primitive blocks and projections](../science/resources-and-regimes.md#2-primitive-blocks-and-success-projections).

### M09. From an overlap to a measured scalar

**Phase:** Background/preliminaries.

**Central question:** How can a coherent matrix element become a numerical answer with controlled failure probability?

**Understanding target:** Explain the Hadamard-test probability, the role of amplitude estimation and the cost of repeating preparations and their inverses. Distinguish the precision of a dimensionless probability from the precision of an energy after rescaling.

**Dependencies and story role:** Uses M01, M07 and M08 to connect coherent access to the classical output contract.

**Conceptual connection:** The scalar-estimation primitive inherits a favorable precision dependence only after its entire controlled preparation has been priced.

**Sources:** Lin §§3.1 and 4.2, especially Example 4.5; [General scalar inverse](../science/scalar-estimation.md#a-general-positive-gap-inverse); [Confidence and output budget](../science/resources-and-regimes.md#5-aggregate-precision-and-confidence).

### M10. Polynomial inverses and normalization

**Phase:** Background/preliminaries.

**Central question:** What must be true before a matrix polynomial is an efficient inverse implementation?

**Understanding target:** Explain how a spectral interval supports an inverse polynomial, how degree controls approximation, and how LCU weights affect circuit normalization. Recognize QSVT and a convergent Neumann polynomial as different implementations with their own preparation and compilation costs.

**Dependencies and story role:** Uses M04 and M07–M09 to connect a matrix interface to inverse-based scalar computation.

**Conceptual connection:** Polynomial approximation becomes useful only when the transformed block, its normalization and its errors fit the same output budget. The project does not equate a physical contraction with a circuit coefficient sum.

**Sources:** Lin §§7.3, 7.5, 8.3 and 8.4; [General inverse](../science/scalar-estimation.md#a-general-positive-gap-inverse); [Finite Neumann implementation](../science/scalar-estimation.md#a-direct-finite-neumann-implementation).

### M11. Classical scalar certificates

**Phase:** Background/preliminaries.

**Central question:** Can a classical trial certify an energy without accurately reconstructing every amplitude?

**Understanding target:** Interpret the positive-definite variational functional and the residual norm bound. Explain why a poor-looking amplitude approximation can still provide an adequate scalar, and why the residual must include every output sector of the full target matrix.

**Dependencies and story role:** Uses M01, M03 and M04 to establish a strong classical alternative before judging the quantum construction.

**Conceptual connection:** The same positive gap that stabilizes a quantum inverse also certifies a classical trial. Output-sensitive comparison applies to both computers.

**Sources:** Lin §4.3 for the linear-system setting; [Classical variational identity](../science/classical-comparison.md#2-positive-definite-variational-error-identity). The scalar functional and residual analysis are supplied by the project, not attributed to the anchor.

### M12. A complete and matched cost question

**Phase:** Background/preliminaries.

**Central question:** What evidence would show that an improved quantum subroutine changes the complete task?

**Understanding target:** Assemble input acquisition, reference/model construction, coherent execution, precision, confidence and scalar output into one ledger. Specify a classical comparison with the same input, target and error. Distinguish a normalization improvement from a proven practical crossover.

**Dependencies and story role:** Integrates M01–M11. This is the global background framework used to judge every later component.

**Conceptual connection:** Scientific computing joins mathematical approximation, data access and numerical certification. A large SD dimension and a fast abstract inverse are insufficient evidence of advantage.

**Sources:** Lin §§1.10, 4.3 and 6.1; [Overview](../science/overview.md); [Classical alternatives](../references/related-work.md#classical-alternatives-that-set-the-comparison); [Resource-account status](../science/resources-and-regimes.md).

## Phase II: narrative and innovation

### M13. The exact SD reduction

**Phase:** Narrative/innovation.

**Central question:** Why is this commutator linear system exactly a Hermitian projected Hamiltonian?

**Understanding target:** Follow how stationarity removes singles from the reference residual and how an excitation acting on a double leaves the SD space or vanishes. Explain the resulting Hermiticity identity and why the same argument cannot simply be extended to unrestricted triples.

**Dependencies and story role:** Uses M02–M04 and M12. This fixes the exact mathematical object before compiling it.

**Conceptual connection:** The Hermitian linear-system interface from Lin is justified by a model-specific excitation-rank argument, rather than imposed on a generally non-Hermitian coupled-cluster problem.

**Sources:** Lin §§6.4 and 8.4 for the required interface; [SD Hermiticity proof](../science/model-and-task.md#why-sd-is-hermitian-here); [Existing LCC framework](../references/related-work.md#existing-quantum-task-and-closest-compiler).

### M14. Neutrality as a projector compression

**Phase:** Narrative/innovation.

**Central question:** Where does the primitive Coulomb norm's volume dependence disappear?

**Understanding target:** Trace the transition synthesis through a discrete divergence, the inverse Laplacian and a second gradient. Explain why their product is an orthogonal projector and yields a volume-independent bound for neutral exchange. Identify this constant-norm mechanism as inherited prior work, with an explicit procedural realization here.

**Dependencies and story role:** Uses M02, M06–M08 and M13. This is the first local compiler mechanism.

**Conceptual connection:** Physical charge cancellation changes the operator being encoded. Fourier access supplies a bounded directional multiplier without requiring globally uniform orbital coefficients.

**Sources:** Lin §§3.3, 6.2 and 6.4; [Neutral projector compression](../science/coulomb-access.md#neutral-exchange-as-a-projector-compression); [Closest compiler and attribution](../references/related-work.md#existing-quantum-task-and-closest-compiler).

### M15. Green columns from a known envelope

**Phase:** Narrative/innovation.

**Central question:** How can a Green-half column be prepared without first computing its global norm?

**Understanding target:** Explain dyadic frequency shells, a known envelope and a bounded amplitude correction. Distinguish the mathematical envelope from the padded uniform-bit circuit envelope, and recognize that the unmeasured failure branches carry the unused norm.

**Dependencies and story role:** Uses M05–M08 and M14 to supply the charged-access primitive that neutral exchange alone does not provide.

**Conceptual connection:** Structured state preparation replaces a large frequency table and an unknown normalization sum with reversible local arithmetic. It does not supply free normalized access to arbitrary data.

**Sources:** Lin §§1.7, 1.8, 2.3 and 6.2 for arithmetic, amplitude and block foundations; [Known-envelope construction](../science/coulomb-access.md#green-half-columns-from-a-known-envelope); [Padded circuit construction](../science/resources-and-regimes.md#2-primitive-blocks-and-success-projections).

### M16. Charged overlaps and their operator

**Phase:** Narrative/innovation.

**Central question:** Which operator is implemented when the two particle positions remain as controls?

**Understanding target:** Distinguish the bounded diagonal two-particle potential from the full Coulomb convolution matrix. Explain why orbital-density columns must be combined coherently with the selector unprepared, and why retaining distinguishable route garbage changes the intended overlap.

**Dependencies and story role:** Uses M08 and M15 to turn the column primitive into the correct direct-interaction block.

**Conceptual connection:** A matrix element and a diagonal operator made from those elements have different norms and input/output registers. Correct operator identification prevents a false constant-normalization claim for the full Coulomb matrix.

**Sources:** Lin §6.2 and §7.3; [Operator encoded by the charged overlap](../science/coulomb-access.md#the-operator-encoded-by-the-charged-overlap); [Primitive blocks and coherent charge mixtures](../science/resources-and-regimes.md#2-primitive-blocks-and-success-projections).

### M17. The single-to-double incidence bound

**Phase:** Narrative/innovation.

**Central question:** Why can a long column still lead to a bounded rank-changing operator?

**Understanding target:** Explain how a source single adds a local pair, why a target double has at most two predecessor singles, and how a bounded squared column norm combines with that predecessor count. Distinguish this argument from flat sparsity counting or a pointwise tail bound.

**Dependencies and story role:** Uses M03, M14–M16 to address the interface between excitation ranks.

**Conceptual connection:** Geometry of the allowed determinant transitions and a Green-kernel norm estimate jointly control the operator; neither ingredient alone completes the construction.

**Sources:** Lin §§6.2–6.3 for access-dependent encoding comparison; [Predecessor proof](../science/singles-doubles-compiler.md#why-a-target-double-has-at-most-two-predecessors); [Uniform source-column bound](../science/singles-doubles-compiler.md#uniform-source-column-bound). The incidence proof is project-specific.

### M18. A coherent compiler for all SD sectors

**Phase:** Narrative/innovation.

**Central question:** How do the separate bounds become one implementable SD block encoding?

**Understanding target:** Explain the row/column maps for the rank-changing block, their retained labels and mutually orthogonal failure sectors. Combine the Fock, direct, exchange and Hermitian off-diagonal blocks with their actual multiplicities. Explain why the dilation need not double the rank-changing normalization.

**Dependencies and story role:** Integrates M13–M17 with M07–M08. This completes the operator interface before preparing its source.

**Conceptual connection:** A norm theorem is converted into a coherent access construction with explicit success conditions, signs and route costs.

**Sources:** Lin §§6.2, 6.4 and 7.3; [Constructive row and column maps](../science/singles-doubles-compiler.md#constructive-column-and-row-maps); [Combining sectors](../science/singles-doubles-compiler.md#combining-the-sectors).

### M19. A source that stays unnormalized

**Phase:** Narrative/innovation.

**Central question:** How does the neutral Gram block prepare the physical source without a hidden success-probability cost?

**Understanding target:** Follow the maximally correlated cell labels, selected exchange block, diagonal mask and coherent ordering step. Explain why the result is the source divided by a known envelope, with failure retained, rather than a normalized source oracle requiring postselection.

**Dependencies and story role:** Uses M08, M14 and M18 to complete the two inputs to the scalar sandwich.

**Conceptual connection:** The same physical structure supplies both operator access and source preparation. A lower bound on the source norm is unnecessary for this construction, although it will matter for demonstrating a nontrivial target accuracy.

**Sources:** Lin §6.2 and §3.1; [Unnormalized source preparation](../science/scalar-estimation.md#an-unnormalized-source-preparation); [Source norm bound](../science/gap-and-signal.md#source-norm-without-a-volume-dependent-preparation-probability).

### M20. A gap certified from the supplied model

**Phase:** Narrative/innovation.

**Central question:** How can the actual SD inverse gap be proved without solving the full matrix?

**Understanding target:** Combine Fock lower bounds with the interaction decomposition. For constant on-site separation, explain the rank-dependent diagonal scaling and the resulting stronger uniform gap. Recognize that the same improved certificate is available to the classical residual calculation.

**Dependencies and story role:** Uses M04, M13 and M18 to replace an abstract inverse promise with a sufficient model-based certificate.

**Conceptual connection:** Preconditioning uses physical excitation energies to sharpen a global spectral bound. It remains conditional on the supplied on-site and hopping bounds; variable on-site energies require an additional term.

**Sources:** Lin §§4.3 and 8.4 for spectral inverse dependence; [General gap](../science/gap-and-signal.md#a-general-sufficient-gap); [Rank-aware gap](../science/gap-and-signal.md#a-stronger-rank-aware-bound-for-constant-on-site-separation); [Classical use of the same bound](../science/classical-comparison.md#all-volume-gap-from-excitation-rank-preconditioning).

### M21. Accuracy that excludes the zero answer

**Phase:** Narrative/innovation.

**Central question:** Is the requested error small enough that a nonzero calculation is necessary?

**Understanding target:** Explain how Fourier aliasing and a variance bound give a lower source signal for the restricted mixing interval and lattice sizes. Connect this to a lower energy magnitude and an additive tolerance that rules out returning zero. Keep this nontriviality certificate separate from computational hardness.

**Dependencies and story role:** Uses M06, M19 and M20 to make the output contract quantitatively meaningful.

**Conceptual connection:** A valid inverse and a bounded source envelope do not alone specify a demanding estimation problem. The lower-signal proof establishes task accuracy, while classical comparison must establish difficulty independently.

**Sources:** Lin §§3.3 and 4.2 for Fourier and estimation foundations; [Uniform lower signal](../science/gap-and-signal.md#a-uniform-lower-signal-for-a-restricted-example-family); [Numerical regime](../science/gap-and-signal.md#numerical-regime-used-throughout). The alias and variance proof is supplied locally.

### M22. Two bounds for one Neumann inverse

**Phase:** Narrative/innovation.

**Central question:** Why are both a physical contraction and a small encoding normalization needed?

**Understanding target:** Distinguish the physical bound called `q` from the circuit normalization called `beta` in the resource program. Explain how the former controls the Neumann tail and the latter controls the geometric LCU sum. Follow the rank-scaled source and fresh success banks in each finite power.

**Dependencies and story role:** Uses M10 and M18–M21 to choose a concrete inverse construction instead of leaving its implementation abstract.

**Conceptual connection:** Convergence is a mathematical property; an efficient coherent realization also needs controlled normalization and paid coefficient preparation. Lin's separate reciprocal-polynomial scaling parameter is not the project's `beta`.

**Sources:** Lin §§7.3, 7.5 and 8.4; [Finite Neumann implementation](../science/scalar-estimation.md#a-direct-finite-neumann-implementation); [Rank preconditioning and explicit inverse](../science/resources-and-regimes.md#4-rank-preconditioning-and-an-explicit-inverse).

### M23. Finite amplitudes from comparisons

**Phase:** Narrative/innovation.

**Central question:** How does a reversible Boolean predicate implement the required bounded amplitude?

**Understanding target:** Explain why unpreparing the uniform comparison register gives an accepted fraction as an amplitude, and why squared predicates are used for square-root factors. Track the zero mode, clipping, fixed-point cosine error and the extra working bits required at small nonzero frequencies.

**Dependencies and story role:** Uses M05, M08 and M15–M16 to replace amplitude-oracle shorthand with finite arithmetic inside the complete construction.

**Conceptual connection:** Arithmetic approximation and coherent selection must be analyzed together; a good top-left block does not imply that an entire approximate dilation is close to an ideal rotation.

**Sources:** Lin §§1.6–1.8 and 6.2; [Bounded amplitudes without inverse square roots](../science/resources-and-regimes.md#3-bounded-amplitudes-without-inverse-square-roots-or-angle-synthesis).

### M24. One coherent error budget

**Phase:** Narrative/innovation.

**Central question:** How do local arithmetic errors become a bound on the final energy?

**Understanding target:** Separate source, matrix, coefficient, finite-polynomial, phase and output-estimation errors. Explain the use of operator/state norms rather than a sum over every SD entry, structural Hermiticity, and comparison with the finite implemented polynomial. Distinguish fixed access bias from physical faults accumulated during execution.

**Dependencies and story role:** Integrates M04, M08 and M19–M23 into a single accuracy contract.

**Conceptual connection:** Lin's error-composition viewpoint becomes a task-level perturbation analysis. Millions of uses do not justify counting the same static matrix approximation again as if it were fresh independent physical noise.

**Sources:** Lin §1.6 and §§7.5, 8.3–8.4; [Finite-precision perturbations](../science/model-and-task.md#finite-precision-perturbations); [Finite arithmetic and Hermiticity](../science/scalar-estimation.md#finite-arithmetic-and-hermiticity); [Aggregate precision](../science/resources-and-regimes.md#5-aggregate-precision-and-confidence).

### M25. Paid output precision and confidence

**Phase:** Narrative/innovation.

**Central question:** How many complete coherent preparations does the chosen scalar estimator actually use?

**Understanding target:** Translate the normalized scalar block into a Hadamard-test probability and then back into energy units. Explain the amplitude-estimation length, median repetitions, inverse preparations, controlled reflections and classical aggregation. Identify why the output stage remains expensive after volume-dependent normalizations have been removed.

**Dependencies and story role:** Uses M09, M19, M22 and M24 to complete the execution-to-output path.

**Conceptual connection:** The quadratic improvement in precision scaling is a property of a full coherent experiment with explicit confidence, not a free measurement of an encoded matrix element.

**Sources:** Lin §§3.1 and 4.2, including Example 4.5; [Aggregate precision and confidence](../science/resources-and-regimes.md#5-aggregate-precision-and-confidence); [Concrete output counts](../science/resources-and-regimes.md#6-concrete-construction-budgets).

### M26. What the concrete resource ledger establishes

**Phase:** Narrative/innovation.

**Central question:** Which consequential improvement is supported by the explicit implementation, and which costs remain unresolved?

**Understanding target:** Connect the degree-11 inverse, finite precision and confidence choices to the conservative logical budget. Explain the roughly 1.85-fold improvement in a scalar normalization envelope. Distinguish that improvement from a measured gate reduction against optimized QSVT, a synthesized circuit, a lower bound, or a fault-tolerant runtime.

**Dependencies and story role:** Integrates M18–M25 with M12. This completes the quantum side of the matched comparison.

**Conceptual connection:** A constructive upper budget makes an abstract algorithm reviewable, while conservative arithmetic, unsynthesized phases and hardware assumptions remain separate from the proved normalization statement.

**Sources:** Lin §§1.7–1.10, 7.3 and 8.4; [Concrete construction budgets](../science/resources-and-regimes.md#6-concrete-construction-budgets); [What improved](../science/resources-and-regimes.md#7-what-improved-and-the-classical-comparison); [Resource verification](../verification/README.md).

### M27. The classical FFT shortcut

**Phase:** Narrative/innovation.

**Central question:** Does spatially varying orbital mixing destroy the cheap classical MP2 scalar?

**Understanding target:** Explain why the unit-dipole kernel and its pointwise square remain convolutions, even when the exchange matrix is noncirculant. Follow the exact FFT evaluation of the diagonal-MP2 scalar with arbitrary allowed mixing coefficients. Identify the remaining bias relative to the full LCCSD target.

**Dependencies and story role:** Uses M06, M11, M14 and M19 to apply the same structural advantages to a classical output calculation.

**Conceptual connection:** Nonuniformity does not erase every convolution. The relevant baseline uses the scalar structure instead of enumerating the full SD space or relying on raw Monte Carlo variance.

**Sources:** Lin §§3.3 and 4.4 for Fourier foundations; [Exact MP2 FFT](../science/classical-comparison.md#3-local-pair-trial-and-an-exact-mp2-fft); [Sampling controls and their limits](../science/classical-comparison.md#7-sampling-controls-when-relevant). The classical derivation is supplied by the project.

### M28. The full residual in a pair representation

**Phase:** Narrative/innovation.

**Central question:** How can a small-support trial certify the unrestricted SD scalar?

**Understanding target:** Identify all three residual output classes: local-pair doubles, singles from the rank-changing block and charge-transfer doubles from hopping. Explain why histories reaching the same determinant must be added before squaring. Connect pair matrices and columnwise FFTs to the complete residual certificate and its work/storage costs.

**Dependencies and story role:** Uses M11, M17–M18, M20 and M27 to turn the cheap trial into a same-target classical calculation.

**Conceptual connection:** Trial support and residual support are different. The classical method retains the action outside its trial space, which is what makes the scalar certificate valid for the full matrix.

**Sources:** Lin §4.3 for the linear-system context; [Complete residual](../science/classical-comparison.md#4-complete-residual-without-the-sd-vector); [Complete classical work](../science/classical-comparison.md#5-complete-classical-work-and-physical-gap-certificates).

### M29. Reading the matched numerical evidence

**Phase:** Narrative/innovation.

**Central question:** What do the finite checks and timings establish about the worked regime?

**Understanding target:** Read the matched energy tolerance, uniform-gap residual bounds and measured times together. Distinguish independent tiny full-fermion checks from larger pair calculations, mathematical identities from their double-precision evaluation, and the size of one pair matrix from peak memory. Recognize that the compact bit-reversal fixture admits smaller-period approximations at reduced precision.

**Dependencies and story role:** Integrates M21, M26 and M28 to assess the actual evidence rather than nominal basis dimension.

**Conceptual connection:** Reproducible calculations test implementation and reveal a strong comparator; they do not prove interval-certified rounding, material accuracy, arbitrary-input hardness or a hardware quantum runtime.

**Sources:** Lin §§1.8 and 1.10 for numerical representation and complexity; [Classical result and interpretation](../science/classical-comparison.md#result-and-interpretation); [Matched input](../science/classical-comparison.md#1-matched-model-and-obtainable-input); [Independent checks](../science/classical-comparison.md#6-independent-structural-and-solver-checks); [Verification limits](../verification/README.md#limits-of-verification).

### M30. The contribution and the remaining advantage gap

**Phase:** Narrative/innovation.

**Central question:** What survives the strongest countercheck, and what would a stronger advantage claim still need?

**Understanding target:** State the result with its qualifiers together: an explicit restricted-model SD compiler with volume-independent repeated-computation normalizations, a concrete inverse refinement, and a matched classical solution in the tested regime. Identify inherited mechanisms, specific implementation work and the incomplete current-version literature comparison. Explain the need for useful obtainable inputs, a demanding scalar regime, strong classical comparisons and compiled quantum resources before asserting advantage.

**Dependencies and story role:** Integrates M12–M29 into the full scientific conclusion. It closes the map by relating the proof, implementation, cost and counterexample evidence.

**Conceptual connection:** A consequential compiler improvement remains scientifically useful when the end-to-end comparison narrows its application claim. The correct endpoint combines the supported construction with the unresolved gap, without promoting a favorable subroutine bound into a speedup theorem.

**Sources:** Lin §§1.10, 4.3 and 8.4; [Overview and scope](../science/overview.md#scope-of-the-result); [Resource comparison and remaining gap](../science/resources-and-regimes.md#7-what-improved-and-the-classical-comparison); [Contribution boundary](../references/related-work.md#what-remains-specific-to-this-construction); [Current-version coverage limit](../references/related-work.md#current-version-coverage-limit).

