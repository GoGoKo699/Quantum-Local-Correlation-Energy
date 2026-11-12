# Matched classical scalar calculation

[Overview](overview.md) · [Quantum resources](resources-and-regimes.md) · [Verification](../verification/README.md) · [Home](../README.md)

## Result and interpretation

The explicit weak-coupling example does not establish quantum advantage. A classical calculation can meet the same requested finite-model energy tolerance without constructing the full singles–doubles amplitude vector.

The strongest implemented comparator is a diagonal-MP2 trial used in the positive-definite variational functional. Its **complete SD residual** can be evaluated with a cell-pair matrix, Fourier transforms, and local hopping. The work is

```math
O(V^2\log V),\qquad \text{storage }O(V^2),
```

for the constant-band-gap example, despite the unrestricted SD dimension being $`\Theta(V^4)`$. The residual includes direct and neutral-exchange terms, the singles produced by $`K^\dagger`$, and the charge-transfer doubles produced by electron and hole hopping. This is a same-target error certificate, not an uncontrolled replacement of LCCSD by MP2.

At the tolerance specified in [the signal certificate](gap-and-signal.md)

```math
\epsilon=5.6980890859\times10^{-9},
```

the evaluated analytical error certificate is below tolerance for all three reported pair-trial instances:

| Cells $`V`$ | Unrestricted SD dimension | Variational $`-e`$ | Energy-error bound, uniform rank gap | Energy-error bound, paid finite gap | Measured total time |
|---:|---:|---:|---:|---:|---:|
| 32 | 247,040 | $`2.4809945\times10^{-6}`$ | $`5.0285262\times10^{-9}`$ | $`4.3791036\times10^{-9}`$ | 0.0006852 s |
| 256 | 1,065,435,136 | $`2.2380260\times10^{-6}`$ | $`4.9025878\times10^{-9}`$ | $`4.4724488\times10^{-9}`$ | 0.0069128 s |
| 2,048 | 4,393,756,786,688 | $`2.2376323\times10^{-6}`$ | $`4.8800545\times10^{-9}`$ | $`4.4947450\times10^{-9}`$ | 0.81534 s |

The times are one run on this execution runtime, including coefficient generation, the periodic Green transform, pair construction, energy, and residual calculation. They are not portable hardware predictions. The certificates use the exact mathematical residual identity below, evaluated in double precision. The implementation is **not an interval-arithmetic certification** of floating-point rounding; that numerical qualification applies to the displayed bounds. The requested tolerance has substantially more margin than the observed small-instance algebraic roundoff, but the latter is not substituted for a formal interval proof.

These examples establish a strong classical comparator and exclude any practical quantum-win claim based on their large unrestricted amplitude dimension. They do not prove that every model admitted by the compiler theorem is classically easy at arbitrary precision. The quantum construction's compact coherent access and normalization result remains a separate statement.

## 1. Matched model and obtainable input

Use the same periodic primitive lattice of side $`L=2^n`$, paired vertically into $`V=L^3/2`$ electronic cells. The primitive interaction is

```math
C=4\pi g(-\Delta_{\rm lat})^+,
```

with zero mode removed, fixed spatial resolution, and the same two-site occupied/virtual orbitals as in the model account. No physical Coulomb-tail replacement is made.

The benchmark has a deterministic compact coefficient program. Index cells lexicographically by $`x\in\{0,\ldots,V-1\}`$, and write $`w=\log_2V`$. Define

```math
p_x=\frac14+
\frac{\operatorname{rev}_w((5x+3)\bmod V)}{2V},
\qquad \gamma_x=\sqrt{p_x(1-p_x)}.
```

Here $`\operatorname{rev}_w`$ reverses the $`w`$ binary digits. The program contains a constant multiplication, addition, reduction, bit reversal, and scaling. It needs no independently supplied disorder table. Since 5 is odd and $`V`$ is a power of two, the coordinate values are a permutation of a dyadic grid in $`[1/4,3/4)`$. The coefficients are spatially nonuniform. This is an explicit algorithm test family, not a sampled material or an assumed hard pseudorandom function. At fixed coefficient precision, discarding sufficiently small output bits permits a smaller-period approximation; the full-precision permutation does not establish intrinsically hard disorder or asymptotic resistance to classical methods.

The remaining parameters are

```math
g=0.01,\qquad E_0=2,
\qquad f^o_{xx}=-E_0/2,\quad f^v_{xx}=E_0/2.
```

For each distinct nearest neighbor, occupied hopping is $`+0.02`$ and virtual hopping is $`-0.02`$, with real symmetric coefficients. Repeated neighbors on very small periodic axes are counted once. Nonzero hopping permits charge transfer and does not preserve the subspace of rigid local excitons. The occupied–virtual Fock block is zero, as required by the specified stationary-reference model.

Every classical computation receives the same coordinate program and electrostatic definition as the quantum algorithm. It explicitly pays to evaluate all coordinates it uses. Constructing a real-space Green kernel requires a finite Fourier transform; a real-space value oracle is not assumed free.

## 2. Positive-definite variational error identity

Write the desired positive scalar as

```math
F=-e=\frac1V b^\dagger A^{-1}b,
\qquad A\succeq\Delta I>0.
```

For any trial amplitude $`t`$, define

```math
\mathcal F(t)=\frac1V
\left(2\operatorname{Re}b^\dagger t-t^\dagger At\right),
\qquad r=b-At.
```

Completing the square gives the exact identity

```math
F-\mathcal F(t)
=\frac1V r^\dagger A^{-1}r.
```

Consequently

```math
0\le F-\mathcal F(t)
\le \frac{\|r\|^2}{V\Delta}.
\tag{1}
```

This is a standard positive-definite variational identity, not a new principle. Its useful application here is that both the trial scalar and the **full residual norm** admit a much smaller representation than the unrestricted SD solution. For the negative correlation energy, (1) gives the interval

```math
e\in[-\mathcal F(t)-\eta,-\mathcal F(t)],
\qquad \eta=\|r\|^2/(V\Delta).
```

The residual certificate controls the same complete finite-model SD scalar. It is not merely a bound on the trial's MP2 value, a local subsystem, or a truncated interaction matrix.

## 3. Local-pair trial and an exact MP2 FFT

Let $`d_x`$ be the unit signed difference of the two primitive sites in cell $`x`$, and define

```math
B_{xy}=d_x^\dagger C d_y,
\qquad J_{xy}=\gamma_x B_{xy}\gamma_y.
```

The matrix $`B`$ is a cell-lattice convolution, while $`J`$ generally is not translation invariant. The source consists of doubles containing two local excitations, with source amplitude $`b_{xy}=J_{xy}`$ for $`x<y`$. The diagonal onsite denominator is $`2E_0`$ for these doubles. Set

```math
T_{xy}=\frac{J_{xy}}{2E_0}\quad(x\ne y),
\qquad T_{xx}=0.
\tag{2}
```

The physical trial amplitude has value $`T_{xy}`$ on the unordered local-pair determinant $`\{x,y\}`$ and zero elsewhere. It has no singles or charge-transfer doubles. Such outputs do nevertheless appear in its residual and must be included.

Even with arbitrary $`\gamma_x`$, the diagonal-MP2 scalar alone is particularly cheap. With $`q_x=\gamma_x^2`$ and $`S_{xy}=B_{xy}^2`$ for $`x\ne y`$, $`S_{xx}=0`$,

```math
F_0=\frac1V\sum_{x<y}\frac{J_{xy}^2}{2E_0}
=\frac{q^\mathsf T S q}{4E_0 V}.
\tag{3}
```

The squared kernel $`S`$ is still a convolution. One FFT convolution therefore computes this scalar in $`O(V\log V)`$ work and $`O(V)`$ storage after or together with the charged Green-kernel construction. Nonuniform orbitals do not remove this classical shortcut. This observation is why a comparison solely against Monte Carlo with raw variance would be weak.

Equation (3) is not generally the full target. The residual calculation supplies the needed correction and bias control.

## 4. Complete residual without the SD vector

All formulas below use real orbitals and the real symmetric hopping matrices of the implemented model. Determinant signs are fixed by canonical increasing electron and hole labels. Signs that multiply one complete output amplitude do not affect its squared norm.

Put $`\kappa_x=2p_x-1`$, and let

```math
U_x=(\rho_{v_xv_x}|C|\rho_{o_xo_x}).
```

In the present two-site model,

```math
U_x=2\gamma_x^2 C_{00}
+(1-2\gamma_x^2)C_{01},
```

where $`C_{00}`$ is the primitive on-site kernel and $`C_{01}`$ the kernel between the paired primitive sites. The direct interaction on a determinant with local excitations at $`x,y`$ is

```math
\mathcal D_{xy}
=\kappa_x\kappa_y B_{xy}-U_x-U_y.
```

Define the symmetric matrix, with its diagonal discarded,

```math
\mathcal I=\mathcal D\odot T+JT+TJ.
\tag{4}
```

The residual separates into three orthogonal output classes.

### Two local excitation pairs

Direct terms preserve the local-pair determinant, while neutral exchange transfers either local pair. Their amplitudes add before squaring. Since the diagonal onsite term cancels $`b`$, the residual for $`x<y`$ is

```math
r^{\rm pair}_{xy}=-\mathcal I_{xy}.
```

Its norm contribution is $`\frac12\sum_{x\ne y}|\mathcal I_{xy}|^2`$.

### Singles from the rank-changing block

The adjoint $`K^\dagger`$ removes either local pair from the trial. Only a local single at $`x`$ survives, with

```math
r^{\rm single}_x
=-\kappa_x\sum_{y\ne x}\gamma_y B_{yx}T_{yx}.
\tag{5}
```

All inserted-pair contributions to the same single are summed coherently. Omitting this channel would give an incorrect certificate even though its numerical norm is small in the reported regime.

### Charge-transfer doubles from hopping

A charged target consists of one local pair at $`y`$ and a separated electron at $`a`$, hole at $`i`$, with $`a\ne i`$, $`y\notin\{a,i\}`$. That target has exactly one local pair, so the three labels identify it uniquely.

There are two histories: move the electron from a trial pair at $`i`$, or move the hole from a trial pair at $`a`$. The complete amplitude, up to one canonical determinant sign, is

```math
f^v_{ai}T_{iy}-f^o_{ai}T_{ay}.
\tag{6}
```

These histories are added before squaring. For the benchmark hoppings the amplitude is $`-h(T_{iy}+T_{ay})`$, $`h=0.02`$. Pauli-forbidden $`y=a,i`$ are removed explicitly.

### Exhaustive norm and scalar

The normal-ordered SD decomposition has no other route from the trial: direct and exchange terms give (4), $`K^\dagger`$ gives (5), and local Fock hopping gives (6). Terms changing excitation rank by two leave the SD space. The three output classes are orthogonal. Hence

```math
\|r\|^2=
\frac12\sum_{x\ne y}|\mathcal I_{xy}|^2
+\sum_x|r^{\rm single}_x|^2
+\sum_{a\ne i}\sum_{y\notin\{a,i\}}
|f^v_{ai}T_{iy}-f^o_{ai}T_{ay}|^2.
\tag{7}
```

The last sum is restricted to nonzero local hopping edges. The variational scalar is

```math
\mathcal F(t)=F_0-
\frac{1}{2V}\sum_{x\ne y}T_{xy}\mathcal I_{xy}.
\tag{8}
```

Equations (7) and (8) are the implemented comparison. They retain every cross-rank and charge-transfer effect needed to certify the full scalar.

## 5. Complete classical work and physical gap certificates

Generate $`p_x`$, the finite Green kernel, and the unit-dipole kernel $`B(r)`$. Then form the pair matrix $`T`$, paying $`O(V^2)`$ arithmetic and storage. Multiplication by $`J`$ is

```math
JX=D_\gamma\,\mathcal F_c^{-1}
\operatorname{diag}(\widehat B)\,
\mathcal F_c\,D_\gamma X,
```

where $`\mathcal F_c`$ transforms the cell coordinate and leaves the column index unchanged. There are $`V`$ columns, giving $`O(V^2\log V)`$ work for $`JT`$. Direct terms, contractions, and hopping residuals cost $`O(V^2)`$ for fixed local degree. The output contains just the scalar and its residual bound.

The implementation retains several pair matrices, coordinate arrays and complex FFT workspaces. The JSON's `pair_matrix_float64_bytes` is the size of **one** $`V\times V`$ matrix, not the peak memory. For $`V=2048`$ one matrix is 32 MiB; the entire calculation uses a larger constant number of such arrays. It does not construct the $`4.39\times10^{12}`$-component SD vector.

### All-volume gap from excitation-rank preconditioning

An all-volume gap suffices for the result. With constant onsite separation $`E_0`$, put $`D=E_0`$ on singles and $`D=2E_0`$ on doubles. For

```math
H=D^{-1/2}(A-D)D^{-1/2},
```

the rank-dependent counts give

```math
\|H\|\le q=
\frac{6(t_o+t_v)+3g\bar G+2\alpha_J+C_K/\sqrt2}{E_0},
```

```math
\bar G=3\pi+\frac{\pi^3}{24},\qquad
\alpha_J=4\pi g\gamma_*^2,\qquad
C_K=4g\gamma_*\sqrt{2\pi\bar G},\qquad \gamma_*=1/2.
\tag{9}
```

Here $`t_o=t_v=h`$ bound the hopping coefficients, and degree six is a common bound for every reported volume. At excitation rank $`m`$, hopping has norm at most $`6m(t_o+t_v)`$ and its denominator is $`mE_0`$. The direct interaction has one primitive charge pair at $`m=1`$ and six at $`m=2`$, each bounded by $`g\bar G`$, so its normalized contribution is at most $`3g\bar G/E_0`$. At most $`m^2`$ exchange lifts contribute, giving at most $`2\alpha_J/E_0`$. The single–double block acquires the factor $`1/(\sqrt2 E_0)`$; its Hermitian off-diagonal dilation has norm $`C_K/(\sqrt2 E_0)`$. Adding these four bounds proves (9).

Since $`q<1`$,

```math
A\succeq(1-q)D\succeq E_0(1-q)I.
```

For the stated parameters, $`q=0.3701902428`$ and
$`\Delta_{\rm rank}=1.2596195145`$. Substituting this all-volume gap in (1) gives the table's uniform-rank certificates, all below the target tolerance. This gap does not use a computed Green diagonal. The rank-preconditioning step requires the constant onsite separation used here; variable onsite energies would need an additional normalized diagonal-deviation bound.

### Optional paid finite gap

The finite Green diagonal is also already available from the paid transform:

```math
G_0(L)=\frac{4\pi}{L^3}\sum_{k\ne0}\lambda(k)^{-1}.
```

The model proof bounds each direct potential by $`gG_0`$ and the rank-changing block using the same diagonal. Therefore the exact finite number can replace the looser all-volume upper bound in that proof:

```math
\Gamma_L=g\left[
6G_0(L)+4\pi+2\sqrt{2\pi G_0(L)}\right],
\qquad
\Delta_L=E_0-2dh-\Gamma_L.
\tag{10}
```

Here $`d`$ is the actual distinct-neighbor degree and $`\gamma_*\le1/2`$. At $`L=8,16`$, the evaluated gaps are 1.380764 and 1.367600. This improvement is paid by the already performed classical input preparation; an equally preprocessed quantum calculation may use it too. It is not a free extra promise supplied only to the classical algorithm.

The saved data retain three bounds with separate names: the older coarse all-volume bound (`Delta_uniform`), the improved all-volume rank bound (`Delta_rank_uniform`), and (10) (`Delta_paid_finite_Green`). The first gives errors near $`7.4\times10^{-9}`$ and does not meet tolerance. Both improved bounds meet tolerance, so computed finite-gap sharpening is not essential to the conclusion. The tighter paid finite bound is included to show the benefit of input preparation already required for the classical scalar.

## 6. Independent structural and solver checks

The small-instance implementation also supplies a general structural SD matrix-vector product. It uses sparse exterior-square Fock actions, diagonal direct terms, and a local-pair incidence map $`E`$. The exchange term on doubles is $`E^\dagger(J\otimes I)E`$; the same map implements $`K,K^\dagger`$. This avoids a dense SD matrix, although it stores the unrestricted SD vector.

At $`L=2`$, an independent second-quantized Hamiltonian is built in the complete 70-dimensional fixed-electron-number space. The structural 52-dimensional SD matrix agrees entrywise to $`8.9\times10^{-16}`$; its source agrees to $`8.7\times10^{-19}`$. This check does not use the structural decomposition to construct its reference matrix.

At $`L=4`$, the structural solver has 247,040 components. Preconditioned conjugate gradients requires nine iterations. The resulting scalar is

```math
F=2.48256483468915\times10^{-6}
```

to the reported floating-point precision. The pair-trial scalar and the three-channel residual norm agree with this structural matrix calculation at $`L=2,4`$: local-pair component discrepancies are below $`2.7\times10^{-19}`$, and residual squared-norm discrepancies below $`2.7\times10^{-23}`$.

The structural calculation also checks diagonal-preconditioned Neumann trials. One correction beyond the diagonal trial has a variational residual bound $`1.31\times10^{-11}`$ at $`L=4`$, using even the looser universal gap. The scalar calculation should therefore stop when its output certificate meets tolerance; solving every amplitude to excessive accuracy is not the fair baseline.

These checks validate algebra and small-instance numerics. They are not material benchmarks, proofs of arbitrary-volume floating-point stability, or fault-tolerant quantum runtime estimates.

## 7. Sampling controls, when relevant

For completeness, the benchmark records exact finite-distribution variances for an MP2 importance sampler. Let

```math
S_B=\sum_{r\ne0}B(r)^2,
\qquad P(r)=B(r)^2/S_B,
\qquad c_B=S_B/(4E_0).
```

Sample $`x`$ uniformly and $`r`$ from $`P`$. Then $`Z=c_Bq_xq_{x+r}`$ estimates (3). For this coefficient program,

```math
\bar q=\frac1V\sum_xq_x
=\frac{11}{48}-\frac1{24V^2}
```

is known analytically. The control

```math
H=c_B\bar q(q_x+q_{x+r}-\bar q)
```

has known mean $`c_B\bar q^2`$. Its residual is

```math
Z-H=c_B(q_x-\bar q)(q_{x+r}-\bar q).
```

At $`L=4`$, the raw variance is $`8.99\times10^{-14}`$, while the centered-pair variance is $`2.86\times10^{-16}`$, a factor of about 314 reduction. The kernel weights and any sampling table still have to be obtained and charged. Since the exact MP2 scalar can already be computed by FFT, these sample variances are not presented as the best runtime. They illustrate why a comparison using raw Monte Carlo variance would be misleading. MP2 sampling also remains a full-target comparator only when its bias or a full-target correction is controlled.

## Reproducibility and scope

The following files are in [verification/](../verification/README.md):

- `benchmark_lccsd.py`: structural SD operator, small-system reference solve, Neumann/variational checks and sampling moments.
- `fermion_reference.py`: independent full fixed-electron-number Hamiltonian builder used only for the tiny validation.
- `pair_variational.py`: pair-matrix scalar and full residual calculation, with no unrestricted SD vector in the reported larger instances.
- `results/benchmark_results.json` and `results/pair_variational_results.json`: observed results.

Run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 python verification/benchmark_lccsd.py --L 2 4
OPENBLAS_NUM_THREADS=1 python verification/pair_variational.py --L 4 8 16
```

The observed runtime used Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 on Linux/x86-64, with `OPENBLAS_NUM_THREADS=1`. No random seed or random input table is used. The exact arithmetic identities and norm proofs establish the theoretical certificates; the saved numerical evaluations use ordinary double precision.

The result is a complete classical countercheck for the stated parameter family and tolerance. It removes the need to speculate about a crossover on these instances. More demanding precision, stronger interactions, different input programs, or first-principles material accuracy require their own matched comparison. None is inferred from an oracle mean-estimation lower bound or from the size of the unrestricted SD basis.
