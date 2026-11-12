# Related work and contribution boundary

[Scientific overview](../science/overview.md) · [Classical comparison](../science/classical-comparison.md) · [Home](../README.md)

This repository specializes established quantum linear-algebra methods to the finite-model LCCSD energy density

```math
e=-\frac{b^\dagger A_{SD}^{-1}b}{V}.
```

Its result is a constructive compiler and a reproducible comparison for an explicit input family. It makes no absolute novelty claim and does not demonstrate quantum advantage. The structured [bibliography](bibliography.json) records the primary sources, versions, and access limits behind this account.

The motivating lead was [Algorithms and compilation](https://github.com/GoGoKo699/Quantum-Off-the-Wall/blob/c2fd6333a0d2f19fd109a10b9baaba6b5c138923/science/algorithms-and-compilation.md), pinned to commit `c2fd6333a0d2f19fd109a10b9baaba6b5c138923`. It served as a discovery guide; the technical claims below rely on primary research and the repository's explicit derivations.

## Existing quantum task and closest compiler

**Chen and Chan [1] supply the existing framework and LCC task.** Their supplemental equation (13) defines the excitation-space commutator linear system; Table II and its accompanying text count input access, matrix normalization, preparation, and correlation-energy readout. Their sparsity-based realization is a useful baseline to improve, not a lower bound on optimized quantum or classical algorithms. This repository restricts the target to singles and doubles and replaces generic access assumptions by a particular procedural lattice model.

**Chen, Sun, and Chan [2] already establish the central constant-norm mechanism.** Appendices C and D explain constant physical norms at fixed excitation rank and give localized transition-pair THC constructions with constant primitive normalizations. Their exchange argument retains angular neutral-density cancellation; their direct factor is a diagonal pair potential. They also warn against artificial normalization from summing momentum blocks. These mechanisms are inherited here. The present specialization additionally constructs the single–double interface, unnormalized source, actual inverse-gap certificate, and complete scalar computation for the specified procedural input.

**Related correlation-energy readout is also prior work.** Rodenas-Ruiz, Zhao, and Lee [11] describe block encodings of quantum-chemical Riccati solutions and correlation-energy-density estimation. Their author-hosted abstract was checked; this repository does not claim amplitude block encoding or scalar energy readout as a new general idea.

## Inherited algorithmic ingredients

| Ingredient | Primary source and checked location | Role here |
|---|---|---|
| Singular-value transformation and inverse approximation | Gilyén, Su, Low, and Wiebe [3], §3.7, Theorem 41 | The general positive-gap inverse is standard QSVT. |
| LCU and projected products | [3], §4.3–4.4, Lemmas 52–53 | Products need distinct success ancillas or explicit intermediate projectors. This is essential for the Gram and source constructions. |
| Coherent scalar estimation | Brassard, Høyer, Mosca, and Tapp [4], §4, Theorem 12 | Amplitude estimation supplies the inverse-precision dependence; it does not make preparation or controlled calls free. |
| Reference distributions, comparisons, and row/column maps | Lemieux et al. [5], §§2–3, 5.3–5.7, Appendix A | Power-law preparation and access-dependent block normalizations are established methods. The repository supplies the concrete shell and incidence maps. |
| Hierarchical kernel encodings | Nguyen, Kiani, and Lloyd [6] | Kernel structure can improve normalization without generic dense loading. A shell or hierarchical encoding alone is not a new mechanism. |

The finite Neumann implementation uses the standard identity $`(I-M)^{-1}=\sum_{k\ge0}M^k`$ when $`\|M\|<1`$. Its application-specific work is to bound the physical contraction $`q`$, the potentially larger circuit normalization $`\beta`$, the source envelope, and all finite errors separately. It introduces no new inverse-series method.

Quantum preconditioning and Green-function computation also have direct precedent in Tong, An, Wiebe, and Lin [13]. The resource ledger's elementary arithmetic is based on conventional reversible addition [14], while optimized function-evaluation circuits [15] are a stronger compilation alternative. Rotation-synthesis algorithms such as Ross and Selinger [16] supply a future compilation route. Their existence does not turn the present analytic operation budgets into a synthesized T count or a hardware runtime.

## Classical alternatives that set the comparison

**Local correlation already exploits long-range cancellation.** Schwilk, Usvyat, and Werner [7] organize strong, close, weak, and distant pairs and improve pair approximations by combining slowly decaying residual contributions before truncation. Neutrality, pair domains, and distant-pair treatment must therefore be available to the classical comparison. Their use of “local coupled cluster” should not be confused with the different abbreviation “linearized coupled cluster.”

**Naive tensor contractions are not the strongest exact baseline.** Mardirossian, McClain, and Chan [8] use alternative real-space, momentum-space, and time representations to reduce exact Gaussian-basis LCCD scaling to $`O(N^5)`$ without a locality approximation. This is a distinct model from the lattice fixture here, but it directly rules out treating a conventional dense contraction count as universally optimal.

**Periodic structure offers classical reductions.** Hansen, Aurbakken, and Pedersen [9] develop periodic local correlation with translation-symmetric tensors in block-Toeplitz form and controlled local domains. The same translation and convolution structure must be retained when comparing quantum and classical work. This repository's periodic support/FFT analysis and its inhomogeneous pair-matrix FFT implementation are explicit matched constructions, not claims that generic FMM or FFT eliminates every coupled-cluster amplitude dimension.

**Modern PNO and related local methods are strong practical competitors.** Pham, Modrzejewski, Nagy, and Klimeš [10] compare local coupled-cluster approaches for molecular solids, including PNO, DLPNO, and LNO treatments and improved long-range approximations. Their molecular-solid calculations are evidence for the strength of practical classical methods, not an exact solver theorem for this repository's finite LCCSD model. The repository consequently tests the same scalar target directly rather than converting a nominal SD dimension into a runtime claim.

For the present fixture, a local-pair trial permits evaluation of both the variational scalar and the **complete** SD residual, including single–double coupling and hopping outside the trial support. For a classical trial vector $`t`$ approximating $`A^{-1}b`$, define

```math
\mathcal F(t)=\frac{2\mathrm{Re}(b^\dagger t)-t^\dagger At}{V},
\qquad r=b-At.
```

With $`A\succeq\Delta I`$, the exact identity

```math
0\le-e-\mathcal F(t)
=\frac{r^\dagger A^{-1}r}{V}
\le\frac{\|r\|^2}{V\Delta}
```

certifies the scalar without storing the full SD vector. The [matched comparison](../science/classical-comparison.md) obtains the required tolerance in the tested regime. Its algebraic certificate and floating-point numerical evaluation are distinguished explicitly. The repository therefore reports a compiler improvement together with an inexpensive classical solution at the stated tolerance; it supplies no supported quantum crossover in that regime.

## What remains specific to this construction

The concrete contribution is the combination of an exact dimer–Poisson model, paid spatially varying coefficient programs, a bounded neutral-response map, coherent Green-half columns with known envelopes, a degree-two single–double predecessor construction, and an unnormalized source for scalar readout. The model also supplies analytical gap and nonzero-signal certificates. Those ingredients remove artificial volume factors from the complete quantum computation under the stated input promises.

They do not establish a first-principles material-accuracy theorem, an unrestricted-triples result, a lower bound on classical algorithms, or a practical quantum crossover. The fixture's compact coefficient description is available to both computers. Its small-period approximations at reduced coefficient precision must not be treated as hard disorder. An improvement over one sparse/envelope implementation is a comparison of constructive upper bounds, not proof that all earlier optimized algorithms require that overhead.

## Current-version coverage limit

Tsemo, Sugisaki, Bhattacharjee, and Prasannaa's *Quantum linear solvers for quantum chemistry: prospects of exponential quantum advantage*, arXiv:2607.08220v2 [12], is a close relevant lead. Bibliographic information and earlier-version material were inspected, but its **current v2 full text could not be retrieved** on 11 November 2026. Direct abstract, HTML, PDF, and export-PDF routes were attempted. The author publication page independently identifies the paper.

The current-version technical comparison therefore remains incomplete. Nothing in this repository should be read as evidence that the inaccessible version lacks an overlapping result, or as an absolute priority claim.

## References

1. **Jielun Chen and Garnet Kin-Lic Chan.** “Framework for robust quantum speedups in practical correlated electronic structure and dynamics.” *Physical Review Letters* **137**, 130601 (2026). [DOI: 10.1103/v2ms-wmz1](https://doi.org/10.1103/v2ms-wmz1). [arXiv:2508.15765v2](https://arxiv.org/abs/2508.15765v2), 22 September 2026; [primary PDF](https://arxiv.org/pdf/2508.15765).

2. **Jielun Chen, Jiace Sun, and Garnet Kin-Lic Chan.** “Quasiparticle quantum simulation of materials with the Bethe–Salpeter equation.” [arXiv:2610.02916v1](https://arxiv.org/abs/2610.02916v1) (2026); [primary PDF](https://arxiv.org/pdf/2610.02916), especially Appendices C–D.

3. **András Gilyén, Yuan Su, Guang Hao Low, and Nathan Wiebe.** “Quantum singular value transformation and beyond: exponential improvements for quantum matrix arithmetics.” *STOC 2019*, pp. 193–204. [DOI: 10.1145/3313276.3316366](https://doi.org/10.1145/3313276.3316366). [arXiv:1806.01838v1](https://arxiv.org/abs/1806.01838v1); [full primary preprint](https://arxiv.org/pdf/1806.01838).

4. **Gilles Brassard, Peter Høyer, Michele Mosca, and Alain Tapp.** “Quantum Amplitude Amplification and Estimation.” *Contemporary Mathematics* **305**, 53–74 (2002). [arXiv:quant-ph/0005055v1](https://arxiv.org/abs/quant-ph/0005055v1); [primary PDF](https://arxiv.org/pdf/quant-ph/0005055v1).

5. **Jessica Lemieux, Matteo Lostaglio, Sam Pallister, William Pol, Karthik Seetharam, Sukin Sim, and Burak Şahinoğlu.** “Quantum sampling algorithms for quantum state preparation and matrix block-encoding.” [arXiv:2405.11436v1](https://arxiv.org/abs/2405.11436v1) (2024); [primary full text](https://arxiv.org/html/2405.11436v1).

6. **Quynh T. Nguyen, Bobak T. Kiani, and Seth Lloyd.** “Block-encoding dense and full-rank kernels using hierarchical matrices: applications in quantum numerical linear algebra.” *Quantum* **6**, 876 (2022). [DOI and primary article](https://quantum-journal.org/papers/q-2022-12-13-876/); [arXiv:2201.11329v4](https://arxiv.org/abs/2201.11329v4).

7. **Max Schwilk, Denis Usvyat, and Hans-Joachim Werner.** “Communication: Improved pair approximations in local coupled-cluster methods.” *Journal of Chemical Physics* **142**, 121102 (2015). [DOI: 10.1063/1.4916316](https://doi.org/10.1063/1.4916316); [author-hosted reprint](https://www.itheoc.uni-stuttgart.de/research/werner/research/files/1.4916316_Reprint.pdf).

8. **Narbe Mardirossian, James D. McClain, and Garnet Kin-Lic Chan.** “Lowering of the complexity of quantum chemistry methods by choice of representation.” *Journal of Chemical Physics* **148**, 044106 (2018). [DOI: 10.1063/1.5007779](https://doi.org/10.1063/1.5007779); [author-institution record and full texts](https://authors.library.caltech.edu/records/avd3x-v0s37).

9. **A. S. Hansen, E. Aurbakken, and T. B. Pedersen.** “Smooth potential-energy surfaces in fragmentation-based local correlation methods for periodic systems.” *Molecular Physics* **119**(9), e1896046 (2021). [DOI: 10.1080/00268976.2021.1896046](https://doi.org/10.1080/00268976.2021.1896046); [primary article](https://www.tandfonline.com/doi/full/10.1080/00268976.2021.1896046).

10. **Khanh Ngoc Pham, Marcin Modrzejewski, Péter R. Nagy, and Jiří Klimeš.** “Sub-kJ/mol Precision at Practical Cost for Molecular Solids through Local Coupled-Cluster Approaches.” *Journal of Physical Chemistry A* **130**(43), 8707–8721 (2026). [DOI: 10.1021/acs.jpca.6c05975](https://doi.org/10.1021/acs.jpca.6c05975); [primary full-text archive](https://pmc.ncbi.nlm.nih.gov/articles/PMC13726384/).

11. **Pablo Rodenas-Ruiz, Andrew Zhao, and Joonho Lee.** “Quantum Solvers for Nonlinear Matrix Equations in Quantum Chemistry.” [arXiv:2605.16189](https://arxiv.org/abs/2605.16189) (2026); [author-hosted abstract and bibliographic record](https://pablorodenas.me/papers/rodenasruiz2026quantum/).

12. **Peniel Bertrand Tsemo, Kenji Sugisaki, Ishita Bhattacharjee, and V. S. Prasannaa.** “Quantum linear solvers for quantum chemistry: prospects of exponential quantum advantage.” [arXiv:2607.08220v2](https://arxiv.org/abs/2607.08220v2) (2026; current full text unavailable in this review); [author publication page](https://sites.google.com/view/srinivasaprasannaa/research-and-publications).

13. **Yu Tong, Dong An, Nathan Wiebe, and Lin Lin.** “Fast inversion, preconditioned quantum linear system solvers, and fast evaluation of matrix functions.” *Physical Review A* **104**, 032422 (2021). [DOI: 10.1103/PhysRevA.104.032422](https://doi.org/10.1103/PhysRevA.104.032422); [arXiv:2008.13295v2](https://arxiv.org/abs/2008.13295v2).

14. **Steven A. Cuccaro, Thomas G. Draper, Samuel A. Kutin, and David Petrie Moulton.** “A new quantum ripple-carry addition circuit.” [arXiv:quant-ph/0410184v1](https://arxiv.org/abs/quant-ph/0410184v1) (2004).

15. **Thomas Häner, Martin Roetteler, and Krysta M. Svore.** “Optimizing Quantum Circuits for Arithmetic.” [arXiv:1805.12445v1](https://arxiv.org/abs/1805.12445v1) (2018).

16. **Neil J. Ross and Peter Selinger.** “Optimal ancilla-free Clifford+T approximation of z-rotations.” *Quantum Information and Computation* **16**(11–12), 901–953 (2016). [arXiv:1403.2975v3](https://arxiv.org/abs/1403.2975v3).
