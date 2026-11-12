# Quantum Local Correlation Energy

This repository studies the complete cost of estimating a linearized coupled-cluster singles-and-doubles (LCCSD) correlation energy density. It gives an explicit quantum compiler for a spatially nonuniform discrete-Poisson electronic model, then tests the same scalar task against a strong classical calculation.

## Conclusion

**The established result is an explicit LCCSD scalar-energy compiler whose operator and source normalizations are bounded independently of volume, with paid procedural input access, explicit preparation, a certified inverse gap and scalar readout costs. This work does not establish end-to-end quantum advantage.**

At the tested parameters and accuracy, the matched classical calculation gives a full-target residual bound below tolerance when evaluated in double precision. This numerical evaluation is not an interval-arithmetic certificate, and it does not establish that every admissible model is classically easy.

## Scientific question

Can locality and neutral transition densities remove a size-dependent compilation bottleneck without hiding the cost in data loading, state preparation, inverse conditioning or output?

The target is one scalar,

```math
e=-\frac{1}{V}b^\dagger A^{-1}b,
```

where $`A`$ is the full singles-and-doubles linearized coupled-cluster matrix and $`V`$ is the number of electronic cells. The model includes full periodic lattice electrostatics, nonuniform orbitals and nonzero hopping. It does not restrict the SD vector to rigid local excitons.

## What is established

- Neutral exchange is encoded through a bounded discrete-gradient projector.
- Charged direct potentials use explicitly prepared Green-half columns with a known constant envelope.
- The single–double interface combines a bounded column norm with at most two predecessors per target.
- The source and scalar readout retain their success projectors and do not require postselection.
- Explicit weak-coupling inequalities certify an actual positive SD gap and a nonzero energy scale.
- A rank-preconditioned Neumann implementation provides a finite inverse construction with a transparent precision and logical-resource ledger.

These mechanisms make the repeated coherent computation depend polylogarithmically on volume when the input coefficients have short reversible programs and the physical parameters and accuracy are fixed. This input restriction is essential. Setup remains separately charged: independent coefficient tables, model construction, reference preparation from material data, and physical approximation error retain their own costs.

The existing LCC framework, fixed-rank neutral-Coulomb ideas, block encodings, LCU, amplitude estimation and inverse approximation are credited in [related work](references/related-work.md). The contribution is the explicit construction and matched cost comparison for the stated model, without an absolute priority claim.

## The classical countercheck

The full SD dimension grows as $`V^4`$, but the requested scalar can be much easier. A local-pair trial and its complete residual can be evaluated in $`O(V^2\log V)`$ work and $`O(V^2)`$ storage. The residual includes both the rank-changing singles channel and charge-transfer doubles from hopping.

For the worked input family, this calculation reaches the requested $`5.6981\times10^{-9}`$ energy-density tolerance at $`V=256`$ and $`V=2048`$, without allocating their billion- or trillion-entry SD vectors. The [comparison](science/classical-comparison.md) gives the exact variational identity, measured numerical results and their floating-point limitations. Large amplitude-vector dimension is therefore insufficient evidence of an advantage.

## Read the account

| Page | Purpose |
|---|---|
| [Overview and result](science/overview.md) | The complete proposition, its consequence and its limits |
| [Model and scalar task](science/model-and-task.md) | Obtainable input, SD identity and finite-precision contract |
| [Coulomb access](science/coulomb-access.md) | Neutral projector and controlled charged-column constructions |
| [Singles–doubles compiler](science/singles-doubles-compiler.md) | Rank-preserving terms and the single–double interface |
| [Gap and signal](science/gap-and-signal.md) | Positive inverse gap and nontrivial requested accuracy |
| [Scalar estimation](science/scalar-estimation.md) | Source preparation, inverse approximation and readout |
| [Resources and regimes](science/resources-and-regimes.md) | Full quantum cost, arithmetic and remaining feasibility gap |
| [Classical comparison](science/classical-comparison.md) | Matched FFT, variational and structural-SD alternatives |
| [Related work](references/related-work.md) | Attribution and the precise comparison boundary |
| [Verification](verification/README.md) | Reproducible scripts, observed outputs and test scope |

For machine-assisted reading, see [llms.txt](llms.txt).

## Learning route

The educational anchor is Lin Lin's [*Lecture Notes on Quantum Algorithms for Scientific Computation*, arXiv:2201.08309v1](https://arxiv.org/abs/2201.08309v1), submitted on 20 January 2022. Start with the [source-to-project bridge](learning/source-bridge.md) for selected sections, notation and the steps supplied by this project. The [30-milestone map](learning/milestones.md) connects 12 background milestones to 18 milestones assembling the compiler, its full cost and the matched classical comparison.

The route explains why encoding a matrix, preparing a source and estimating one energy density have different costs. It ends with the supported restricted-model result and the remaining gap to advantage.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

The repository uses the [MIT license](LICENSE).
