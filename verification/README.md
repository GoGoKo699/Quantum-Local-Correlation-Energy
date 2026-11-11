# Reproducing the algebra, resources and classical comparison

[Home](../README.md) · [Scientific overview](../science/overview.md) · [Classical comparison](../science/classical-comparison.md) · [Quantum resources](../science/resources-and-regimes.md)

The scripts check the explicit finite model and evaluate the stated resource formulas. The classical programs perform numerical computation. The quantum ledger evaluates an analytical circuit construction; it does not compile or run a quantum circuit.

## Environment and commands

The recorded numerical runs use Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 on Linux/x86-64. The package versions are recorded in [requirements.txt](requirements.txt). The resource calculator and its rational-arithmetic validation use only the Python standard library.

From the repository root, install the numerical dependencies in your chosen Python environment:

~~~sh
python -m pip install -r verification/requirements.txt
~~~

Then run the following commands. Each script writes its result into the results directory; the resource calculator takes an explicit output path.

~~~sh
OPENBLAS_NUM_THREADS=1 python verification/check_poisson_algebra.py
OPENBLAS_NUM_THREADS=1 python verification/check_signal.py
OPENBLAS_NUM_THREADS=1 python verification/benchmark_lccsd.py --L 2 4
OPENBLAS_NUM_THREADS=1 python verification/pair_variational.py --L 4 8 16
python verification/resource_ledger.py --output verification/results/resource_results.json
python verification/validate_resource_ledger.py
~~~

The structural SD solver is intended for the small sizes shown. It stores the full SD vector, whose length grows as the fourth power of the number of cells. The larger pair-variational instances use a different representation and avoid that vector.

The reported wall times are observations from individual runs, including the setup and computation described in the corresponding script. They are not a portable CPU specification, a statistical timing study, or a comparison in seconds with an unbuilt quantum computer.

## What each check establishes

| Program | Independent object or property checked | Saved output |
|---|---|---|
| [check_poisson_algebra.py](check_poisson_algebra.py) | Full fixed-electron Hamiltonian against the procedural Poisson, neutral Gram, source and single–double formulas; norm and gap inequalities at a tiny size | [poisson_algebra.json](results/poisson_algebra.json) |
| [check_signal.py](check_signal.py) | Primitive Fourier kernel against the aliased cell symbol; frequency fractions and signal bounds at $L=8,16,32$ | [signal_certificate.json](results/signal_certificate.json) |
| [benchmark_lccsd.py](benchmark_lccsd.py) | Structural SD matrix against an independent second-quantized builder at $L=2$; iterative scalar solves, residual bounds, MP2 FFT and controlled expansions | [benchmark_results.json](results/benchmark_results.json) |
| [pair_variational.py](pair_variational.py) | Complete SD residual of the local-pair trial, validated against structural SD at $L=2,4$, then evaluated without an SD vector at $L=4,8,16$ | [pair_variational_results.json](results/pair_variational_results.json) |
| [resource_ledger.py](resource_ledger.py) | Explicit inverse degree, normalization, word lengths, error/confidence allocation, primitive-call counts and conservative gate/width formulas | [resource_results.json](results/resource_results.json) |
| [validate_resource_ledger.py](validate_resource_ledger.py) | Padded-shell identities, squared-predicate amplitude counting, fixed-point cosine values against rational enclosures, and aggregate error budgets | [resource_validation.json](results/resource_validation.json) |

[fermion_reference.py](fermion_reference.py) is the shared tiny-system reference builder. It enumerates the complete fixed-electron sector and is used only to check the formulas independently. It is not the scalable classical comparison.

## Numerical evidence

The independent structural matrix comparison at $L=2$ has maximum entry discrepancy below $9\times10^{-16}$. The source discrepancy is below $9\times10^{-19}$. For the pair method, the complete residual norm agrees with the structural SD calculation at $L=2,4$ to below $3\times10^{-23}$ in squared-norm difference. The primitive/cell Fourier alias agrees within $1.5\times10^{-15}$ at the checked signal sizes.

The uniform rank-aware gap is $1.2596195145$ for the worked parameters. With that gap, the pair method's evaluated full-target energy error bounds are approximately $4.9026\times10^{-9}$ at 256 cells and $4.8801\times10^{-9}$ at 2,048 cells, below the target $5.6981\times10^{-9}$. The finite Green diagonal supplies an optional sharper bound, recorded separately.

The quantum calculation uses degree 11, 26-bit bounded-amplitude comparisons and 24-bit finite-LCU intervals at that target. Its explicit coherent amplitude-estimation construction uses $2^{20}$ as its phase-estimation parameter and eleven independent repetitions, giving an evaluated failure upper bound below one percent. All arithmetic, confidence and call-count assumptions are exposed in the ledger.

## Limits of verification

The uniform analytical arguments are in the scientific pages. Finite numerical agreement supports their implementation; checking several lattice sizes does not prove an asymptotic statement.

The classical residual identity is exact, but its displayed numerical evaluation uses ordinary double precision. The saved residual bounds are not interval-arithmetic guarantees for the complete numerical calculation. Rational enclosures in the quantum primitive checker certify the tested cosine evaluations only, together with the separately stated uniform approximation argument.

No script establishes LCCSD accuracy against an exact material ground state, a classical complexity lower bound, practical quantum advantage, a compiled Clifford+T count, fault-tolerant physical resources, or hosted mathematical rendering.
