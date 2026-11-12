# A paid logical-resource account for procedural Poisson LCCSD

[Overview](overview.md) · [Scalar estimation](scalar-estimation.md) · [Classical comparison](classical-comparison.md) · [Verification](../verification/README.md) · [Home](../README.md)

**Status, 11 November 2026.** This note specifies a complete scalar estimator,
including an explicit inverse polynomial, preparation, cleanup, output,
precision, and confidence. Its executable ledger gives **analytic construction
budgets, not counts from a compiled or simulated quantum circuit**. Rotation
synthesis and fault tolerance have not been run. The matched classical examples
are easy at the stated tolerance; no practical or asymptotic classical
separation is established.

The useful refinement is to encode the rank-preconditioned residual directly.
For the stated instance its physical norm is at most 0.370191 and its encoding
normalization is 0.617265. A degree-11 Neumann polynomial therefore supplies a
fully specified inverse at the chosen accuracy. The distinction between those
two constants, and the cost of coherent amplitude estimation, matters more than
the formal dimension of the singles/doubles matrix.

## 1. Task, input, and what is charged

The finite model, access construction and underlying SD block theorem are in
[the supplied model](model-and-task.md),
[the Coulomb access construction](coulomb-access.md),
[the gap certificate](gap-and-signal.md),
and [the SD compiler](singles-doubles-compiler.md).
Here $L=2^n\ge4$, the primitive torus has $N=L^3$ sites, and there are
$V=L^3/2$ electronic cells. Two sites in each cell support exactly orthonormal
occupied and virtual orbitals, with density weights $p_x,1-p_x$ and neutral
transition density

$$
T_x=\gamma_x(\delta_{s_x+e_z}-\delta_{s_x}),\qquad
\gamma_x=\sqrt{p_x(1-p_x)}\le\gamma_*=1/2.
$$

The complete finite periodic interaction is

$$
C=4\pi g(-\Delta_{\rm lat})^+,\qquad
\lambda(k)=4\sum_{j=1}^3\sin^2(k_j/2),
$$

with zero mode removed. There is no distance cutoff. The target is the single
real number $e=-b^\dagger A_{SD}^{-1}b/V$. Neither the full solution vector nor
all matrix entries are output.

The priced instance uses exactly the same input program as the classical
benchmark:

$$
p_x=\frac14+\frac{\operatorname{bitreverse}_{w}((5x+3)\bmod V)}{2V},
\quad w=\log_2V=3n-1,
$$

with $g=.01,E_0=2,t_o=t_v=.02$. The occupied/virtual onsite entries are
$-E_0/2,+E_0/2$, and the real nearest-neighbor hoppings have magnitudes .02.
The $p_x$ field is an exact $3n$-bit rational, computed by constant arithmetic
and a wire reversal. Arbitrary disorder tables do not inherit this cost.
Nonuniform bounded hoppings could be used by adding their reversible program
cost, but the numerical ledger prices the specified constant hopping program.

The effective model is supplied. Deriving an appropriate Fock operator or local
orbitals from nuclei, and assessing this model against a real material, remain
separate tasks. Small classical coefficient tables, arithmetic compilation,
and phase compilation are setup costs. They are described below and are not
timed. There is no QRAM assumption.

## 2. Primitive blocks and success projections

The neutral map and Gram block are

$$
W=C^{1/2}T
 =\sqrt{4\pi g}\,\operatorname{diag}
 \frac{e^{-ik_z}-1}{\sqrt{\lambda(k)}}\,F S M_\gamma,
\qquad J=W^\dagger W.
$$

Their normalizations are

$$
\alpha_W=\sqrt{4\pi g}\gamma_*,\qquad
\alpha_J=4\pi g\gamma_*^2.
$$

Every expression here denotes a **projected-unitary block**. For example, the
Gram block uses $U_W^\dagger\Pi_W U_W$, implemented with separate signal
registers or an equivalent explicit intermediate success projection. The plain
unitary product $U_W^\dagger U_W$ is the identity and is incorrect. Clean
arithmetic workspace can be reused after its exact inverse; signal registers
that carry failures are retained until the final projection.

From a maximally entangled pair of cell labels, the $J$ block prepares its
vectorized matrix. Masking collisions, coherently sorting $x,y$, and projecting
the order qubit gives

$$
B=\frac{b}{\sqrt V B_0},\qquad B_0=\alpha_J/\sqrt2.
$$

The factor $\sqrt2$ uses a real symmetric Gram matrix. Our construction makes
that property structural: if finite phase synthesis produces a complex
Hermitian Gram matrix, coherently average its circuit and its entrywise
conjugate. The normalization is unchanged. The ledger charges four $W$
invocations for this real Gram source, rather than assuming the two branches are
free. The inverse circuits and conjugate circuits are literal inverses and
conjugates of the chosen finite gate strings.

### Explicit half-Green columns with $Z=8\pi$

This implementation slightly enlarges the earlier $7\pi(1-1/L)$ envelope in
exchange for simple uniform-bit preparation. Let

$$
B_s=[-2^s,2^s)^3,\quad S_0=B_0\setminus\{0\},\quad
S_s=B_s\setminus B_{s-1}\ (s\ge1).
$$

Prepare a uniform $n$-bit integer $q\in[0,L)$. The zero integer is a failure
branch. For $q\ne0$, compute the highest-set-bit index $s$. There are $2^s$
integers in that shell branch. Their lower $s$ uniform bits can be erased by
controlled Hadamards, with the leading bit and bookkeeping retained or cleared
reversibly. Thus shell $s$ has unnormalized probability $2^s/L$.

For $s=0$, three uniform bits select the eight points of $B_0$; the origin
is a failure. For $s\ge1$, six bits select 64 equal cube labels, of which 56
are the outer cubes of side $2^{s-1}$. The eight inner cubes are failures.
Uniform low coordinate bits complete a point in the selected cube. The shell
label is recoverable from every good frequency and is erased there. Invalid
branches retain their success flags.

Every good frequency in $S_s$ has squared reference amplitude

$$
\frac{1}{8L4^s}.
$$

The bounded correction $r_s(k)=2^{s+1}/(L\sqrt{\lambda(k)})\le1$ gives

$$
q_r=\frac1{\sqrt{2L^3}}
 \sum_{k\ne0}\frac{e^{-ik\cdot r}}{\sqrt{\lambda(k)}}|k\rangle,
\qquad g(8\pi)\langle q_r|q_s\rangle=C_{rs}.
$$

The squared success norm is $G_0/(8\pi)$, where
$G_0=(4\pi/N)\sum_{k\ne0}\lambda(k)^{-1}$. Its value is not needed for
preparation. There is no measurement or amplitude amplification of this
success. The overlap must include the success projector between its two
preparations. With the position labels retained, it encodes a diagonal
two-particle potential, not the full single-particle matrix $C$ with constant
normalization.

For a charge density, the two primitive columns must be combined coherently.
Prepare a uniform $3n$-bit integer $j$, select site $s_x$ when
$j<2^{3n}p_x$, select $s_x+e_z$ otherwise, then uncompute the selector and
unprepare $j$. Projecting this uniform register gives
$p_xq_{s_x}+(1-p_x)q_{s_x+e_z}$, with normalization one. Retaining its route
register instead would give the wrong Gram matrix. A uniform sign branch
forms an electron-minus-hole density, with normalization two.

For the single/double block, that difference column is followed by $W^\dagger$
and the canonical occupation/sign masks. At fixed source single, the map from
new local pair label to target double is injective. A target double has at most
two predecessor singles. Their row isometry is combined with the retained
source labels only at the row/column overlap. Thus

$$
\alpha_K=4g\gamma_*\sqrt{2\pi Z},\qquad Z=8\pi.
$$

The Hermitian off-diagonal dilation uses $K$ and the literal adjoint of the
same implemented block. It does not independently round the two directions.

## 3. Bounded amplitudes without inverse square roots or angle synthesis

Let $u=j/2^b$ for a uniform $b$-bit integer $j$. Compute a Boolean predicate
into a retained flag, undo its arithmetic, apply the inverse Hadamards to the
uniform register, and include **both** that register and the flag in the
success projector. The good amplitude is the fraction of accepted integers,
not its square root.

This implements the needed bounded amplitudes with squared predicates:

| Factor | Predicate |
|---|---|
| $\gamma_x/\gamma_*$, for $\gamma_*=1/2$ | $u^2<4p_x(1-p_x)$ |
| $\lvert e^{-ik_z}-1\rvert/\sqrt\lambda$ | $u^2\widehat\lambda<\widehat\lambda_z$ |
| Shell correction | $u^2L^2\widehat\lambda<4^{s+1}$ |

The zero mode is explicitly rejected; approximate numerators are clipped into
their allowed interval. The strict inequality gives a ceiling convention at
nonintegral thresholds, with absolute discretization error at most $2^{-b}$.
It handles zero and unit endpoints. The phase of the neutral multiplier is
separate: $e^{-ik_z}-1=-2i e^{-ik_z/2}\sin(k_z/2)$, so only a sign and dyadic
linear phase are needed.

Evaluate $\lambda=2\sum_j(1-\cos k_j)$ and
$\lambda_z=2(1-\cos k_z)$. If each cosine has error at most $\delta$, the
errors in those quantities are at most $6\delta,2\delta$. On nonzero modes,
$\lambda\ge16/L^2$. Once $\widehat\lambda\ge\lambda/2$, the clipped ratio
has error at most $\delta L^2$. The inequality
$|\sqrt a-\sqrt b|\le\sqrt{|a-b|}$ therefore gives amplitude error at most
$L\sqrt\delta$. The shell ratio obeys the same conservative bound. We choose

$$
\delta=2^{-s},\quad s=2b+2n+10,
$$

so this contribution is at most $2^{-b-5}$. This intentionally spends more
bits than a carefully conditioned direct inverse-square-root circuit.

For a fully specified arithmetic construction, reduce the cosine argument to
$[0,\pi/2]$, with its sign retained. Use a degree-$d$ polynomial in $t=x^2$,
$P_d(t)=\sum_{j=0}^d(-1)^jt^j/(2j)!$. Choose the smallest $d$ with

$$
\frac{1.6^{2d+2}}{(2d+2)!}\le 2^{-(s+2)}.
$$

Horner evaluation stores its previous iterates for later uncomputation.
Coefficients and products are truncated with

$$
f=s+\lceil(d+2)\log_2 3\rceil+8
$$

fractional bits. On the reduced interval, $t<3$; summing the geometric
propagation of coefficient and multiplication errors bounds it by
$(3^{d+1}+8)2^{-f}$. The Taylor remainder and the error of a sufficiently
accurate fixed constant $\pi$ fit inside the remaining $2^{-s}$ budget.
Exact integer products perform the final squared comparisons. No derivative
bound at $\gamma_x=0$ is required, since $p_x$ is exact and the gamma
predicate uses its squared value directly.

An approximate comparator dilation need not be close to a chosen ideal
rotation as a full unitary. The analysis uses the actual success block and its
operator perturbation, followed by a structurally Hermitian matrix. It does
not infer full-unitary closeness from a close top-left entry.

## 4. Rank preconditioning and an explicit inverse

Let $D=E_0 I$ on singles and $D=2E_0 I$ on doubles. Define

$$
M=-D^{-1/2}(A-D)D^{-1/2},\qquad
A=D^{1/2}(I-M)D^{1/2}.
$$

This sign convention makes the inverse polynomial $\sum_k M^k$. Constant
onsite separation is part of this specialization. Variable onsite energies
require adding their scaled residual to $M$ and to its bounds.

We encode $M$ directly. In rank $m=1,2$, the hopping bound is
$m6(t_o+t_v)$, the number of direct pairs is $\binom{2m}{2}$, and the number
of exchange choices is $m^2$. Division by $mE_0$, and the off-diagonal
denominator $\sqrt2 E_0$, gives

$$
\boxed{\beta=
 \frac{6(t_o+t_v)}{E_0}+\frac{3gZ}{E_0}
 +\frac{2\alpha_J}{E_0}+\frac{\alpha_K}{\sqrt2 E_0}.}
$$

This is an **encoding normalization**, not the physical norm. Replacing $Z$
by $\bar G=3\pi+\pi^3/24$ and $\alpha_K$ by
$4g\gamma_*\sqrt{2\pi\bar G}$ gives the independent physical certificate
$\lVert M\rVert\le q$. The example yields

| Quantity | Value |
|---|---:|
| Physical contraction $q$ | 0.370190242753 |
| Encoding normalization $\beta$ | 0.617264703730 |
| Certified gap $E_0(1-q)$ | 1.259619514494 |
| Earlier unpreconditioned gap certificate | 0.827218050007 |

Both $q<1$ and $\beta<1$ matter. Physical convergence alone does not make
this geometric LCU efficient: if $\beta>1$, its coefficient normalization
can grow exponentially in the truncation degree.

Because $b$ lies in doubles, $c=D^{-1/2}b=b/\sqrt{2E_0}$, and
$B_c^2=B_0^2/(2E_0)$. The already prepared normalized-by-envelope source is
unchanged: $c/(\sqrt V B_c)=B$. Define

$$
F_p=\sum_{k=0}^pM^k,\qquad
L_p=\sum_{k=0}^p\beta^k\le(1-\beta)^{-1}.
$$

Its LCU uses weights $\beta^k/L_p$ and a controlled product of $k$ copies of
the $M/\beta$ block. Each copy has **fresh success ancillas**. Repeating a
unitary using the same flags does not encode a matrix power. One safe SELECT
applies the $j$th block when the degree label satisfies $k\ge j$, for
$j=1,\ldots,p$. Arithmetic workspace is reused; all earlier signal banks are
retained. The output scalar sandwich likewise uses distinct source-success
banks for $B$ and $B^\dagger$.

The $p+1$ coefficients are a small, paid classical table. A uniform $h$-bit
integer and cumulative interval comparisons select $k$; the selector and
uniform register are unprepared after SELECT. Dyadic interval weights give
$L_1$ coefficient error at most $2(p+1)2^{-h}$. This charges both comparison
directions and all uniform Hadamards, and requires no free amplitude oracle.
The same technique prepares the constant-size rank and interaction choices.

Preconditioning and quantum inverse/Green-function computation are inherited
mechanisms; see Tong, An, Wiebe and Lin,
[Phys. Rev. A 104, 032422 (2021)](https://arxiv.org/abs/2008.13295).
The rank factors and the bounds above specialize them to this supplied model.
Constant fixed-rank Coulomb normalization is also prior work, explicitly
credited in the underlying notes to
[Chen, Sun and Chan, arXiv:2610.02916](https://arxiv.org/abs/2610.02916).

## 5. Aggregate precision and confidence

The ideal polynomial truncation obeys

$$
|e-e_p|\le B_c^2\frac{q^{p+1}}{1-q}.
$$

Let the implemented Hermitian block be $\widetilde M$, with error $d_M$,
and let the implemented source differ by $d_B$. Both source vectors have
norm at most one as projected preparations. Put $r=q+d_M<1$. Comparing the
**finite** polynomials first, by
$\lVert\widetilde M^k-M^k\rVert\le k r^{k-1}d_M$, and then adding the ideal tail
gives the sufficient implementation bound

$$
E_{\rm impl}\le B_c^2\left[
 \frac{d_M}{(1-r)^2}+\frac{2d_B}{1-r}
 +2L_p(p+1)2^{-h}\right].
$$

There is no omitted tail for a perturbed infinite inverse. This comparison
never substitutes that inverse for the finite implemented polynomial.

The ledger conservatively takes $\eta=2^{-b}$,
$d_M=64\beta\eta,d_B=32\eta$. These constants cover the two amplitude
predicates and phase transform in $W$, projected Gram products, coherent
charge mixtures, K, real averaging, and finite constant-size LCU coefficients.
For example a $W$ block has error below $5\eta$, its Gram below
$10\eta$, a half-Green charge column below $3\eta$, and the real source
below $\sqrt2\,10\eta$; the displayed allowances are larger. Top-level
constant coefficients use $b+8$ interval bits. Dyadic Fourier phase rotations
receive a separate per-primitive operator-error allocation below $\eta$.

Allocate $\epsilon/4$ to the ideal tail, $\epsilon/8$ to primitive error,
$\epsilon/8$ to geometric coefficient error, and $\epsilon/2$ to output
estimation. A Hadamard test of the complete projected scalar unitary has
probability

$$
a_H=\frac{1+\operatorname{Re}\langle\widetilde B|
             \widetilde F_p/L_p|\widetilde B\rangle}{2}.
$$

The reported energy is $-B_c^2L_p(2\widehat a_H-1)$. Therefore it suffices to
estimate $a_H$ to $\epsilon/(4B_c^2L_p)$.

The explicit output algorithm used here is the original coherent amplitude
estimator of Brassard, Høyer, Mosca and Tapp,
[Theorem 12](https://arxiv.org/pdf/quant-ph/0005055). For power-of-two $M$,
one run has probability at least $8/\pi^2$ of error at most

$$
\frac{2\pi\sqrt{a_H(1-a_H)}}{M}+\frac{\pi^2}{M^2}
\le\frac\pi M+\frac{\pi^2}{M^2}.
$$

It uses $M-1$ Grover iterations and $2M-1$ Hadamard-test preparations or
their inverses. Controlled Grover iterations need control only their
reflections: the surrounding preparation and its inverse cancel on the
inactive control branch. Those preparations already include a controlled
scalar circuit for the Hadamard test. The ledger charges controlled
zero-state/good-state reflections and the inverse QFT of the AE counter.

For confidence, take the median of an odd number $R$ of independent scalar
estimates. The calculator evaluates the binomial upper tail using per-run
failure at most $1-8/\pi^2+2\cdot10^{-6}$; the final term allows an AE-counter
inverse-QFT operator error of $10^{-6}$. The example uses $R=11$, giving
failure at most 0.008875. Classical aggregation costs $O(R\log R)$ and only
one final scalar is read out.

## 6. Concrete construction budgets

Run:

```sh
python verification/resource_ledger.py \
  --output verification/results/resource_results.json
python verification/validate_resource_ledger.py
```

At $\epsilon=5.698089085861593\times10^{-9}$, the all-volume result is:

| Item | Explicit choice or evaluated bound |
|---|---:|
| $B_c^2$ | $1.23370055013617\times10^{-4}$ |
| Neumann degree $p$ | 11 |
| $L_p$ | 2.60477792698345 |
| Amplitude bits $b$ | 26 |
| Geometric interval bits $h$ | 24 |
| Ideal truncation error | $1.29748\times10^{-9}$ |
| Primitive implementation error | $3.69900\times10^{-10}$ |
| Geometric coefficient error | $4.59697\times10^{-10}$ |
| Output estimation error | $1.92558\times10^{-9}$ |
| Total additive error bound | $4.05266\times10^{-9}<\epsilon$ |
| AE $M$, repetitions $R$ | $2^{20},11$ |
| Hadamard-test preparations/inverses | 23,068,661 |
| $W$ invocations, padded upper count | 2,214,591,456 |
| Charge half-Green invocations, padded upper count | 2,030,042,168 |

These large output counts are a constructive upper algorithm, not a lower
bound. Better amplitude-estimation constants or scalar estimators may improve
them. The earlier dimensionless factors that omitted logarithms, confidence,
and access were not gate counts and should not be compared directly to these.

### Arithmetic and routing convention

The gate ledger deliberately uses a simple reversible schoolbook construction.
It pads a $w$-bit ripple adder to $9w+8$ elementary Boolean gates. This
follows from the Toffoli/CNOT/NOT counts of
[Cuccaro, Draper, Kutin and Moulton](https://arxiv.org/pdf/quant-ph/0410184).
A compare uses two additions plus $4w+8$ gates. A signed full multiplication,
implemented by sign extension and controlled shifted additions, is bounded
by $64w^2+64w$ such gates. Truncated output copying and multiplication
uncomputation are charged separately. This is intentionally less efficient
than specialized arithmetic such as
[Häner, Roetteler and Svore](https://arxiv.org/abs/1805.12445).

The uniform word width is
$w_{\mathrm{arith}}=f+2b+6n+16$; it covers the exact squared comparisons, including
the $3n$-bit input field. The program exposes every cosine, multiplier,
comparator and LCU term used in its majorant. Shell routing is bounded by
$16(n+1)$ word comparisons and $64(n+1)$ word operations. SD routing is
padded to 128 comparisons, 256 word controlled swaps and 48 word additions.
These cover the constant-sized choices of up to four quasiparticles, six
pairs, up to 24 neighbor moves, collision masks, canonical sorting, inverse
routes, and internal constant-size interval selectors. No routing over $V$
or over the SD dimension is enumerated.

The circuit-level majorant is

$$
C_M\le8C_q+8C_W+C_{\rm rank}+C_{\rm topLCU},\qquad
C_B\le4C_W+C_{\rm rank}+C_{\rm pair},
$$

$$
C_{\rm scalar}\le2C_B+pC_M+C_{\rm finiteLCU}+C_{\rm degree\ controls}.
$$

The direct branch is charged up to four $q$ calls; exchange up to four $W$;
the real Hermitian K dilation up to four of each. These deliberately count
coherent alternatives separately. Six additional coherent controls are padded
onto **every** primitive gate. A controlled conjunction implemented with clean
ancillas gives at most 13 Toffolis per Boolean gate at this control depth;
controlled Hadamards additionally use at most two exact T gates, and controlled
dyadic phase gates at most three single-qubit phases. Reflections receive a
separate width-dependent allowance. This is a conservative construction-level
majorant, not a gate-minimized netlist.

| $L$ | $V$ | Logical qubits, analytic upper budget | Toffolis, analytic upper budget | Unsynthesized phase rotations |
|---:|---:|---:|---:|---:|
| 4 | 32 | 37,966 | $1.96256\times10^{19}$ | $1.19588\times10^{11}$ |
| 8 | 256 | 41,514 | $2.14351\times10^{19}$ | $2.57446\times10^{11}$ |
| 16 | 2,048 | 45,752 | $2.51698\times10^{19}$ | $4.51777\times10^{11}$ |
| 64 | 131,072 | 52,896 | $2.94821\times10^{19}$ | $1.00985\times10^{12}$ |
| 256 | 8,388,608 | 60,946 | $3.69539\times10^{19}$ | $1.79382\times10^{12}$ |

The width includes $p$ fresh $M$-signal banks, two source banks, the degree
and uniform interval registers, shared arithmetic workspace, and the AE
counter. It makes no ancilla-compression claim. The large gate budgets reflect
full-width multiplication, retained Horner history, generic controls and
conservative repetitions; they are not evidence that an optimized circuit
must use this many gates.

Each access circuit's phase synthesis is allocated an operator error per phase
of at most $\eta/(32R_{\rm phase})$, where $R_{\rm phase}$ bounds the number
in a $W$ or $q$ primitive. AE-counter phases have their separate allocation.
An eventual Clifford+T count can use seven T gates per exact Toffoli, the
displayed exact controlled-H overhead, and the **actual synthesized** phase
counts. We do not replace those counts by a guessed constant times a logarithm.
[Ross and Selinger](https://arxiv.org/abs/1403.2975) give applicable synthesis
methods and a typical asymptotic cost, not the finite counts of this uncompiled
circuit. Controlled phase compilation must preserve its specified phase,
including phase corrections relevant when controlled.

Once those finite gate strings are fixed, AE acts on their exact logical
unitaries. Their access error is already a perturbation of the estimated
matrix/source. It need not be multiplied again by the millions of AE calls.
Logical faults, noisy rotations, and physical failure budgets would be an
additional error analysis; no physical runtime or hardware qubit estimate is
provided here.

## 7. What improved, and the classical comparison

The rank-aware construction replaces an unspecified inverse transformation by
a concrete degree-11 inverse polynomial and improves the all-volume gap
certificate from 0.827218 to 1.259620. The infinite geometric scalar envelope
$B_c^2/(1-\beta)\simeq3.22338\times10^{-4}$ is about 1.85 times smaller than
$B_0^2/\Delta\simeq5.9655\times10^{-4}$ from the earlier gap certificate.
This is a modest useful normalization improvement. It is not a demonstrated
gate-count reduction over an optimized, fully compiled QSVT inverse: no such
baseline phase sequence has been compiled here. Standard QSVT remains an
alternative when its polynomial/phase synthesis and robustness are charged;
see [Gilyén, Su, Low and Wiebe](https://arxiv.org/abs/1806.01838).

Classical algorithms receive the identical procedural field and FFT Green
kernel. With constant onsite denominators, the MP2 scalar for arbitrary
$\gamma_x$ reduces to one FFT convolution of $\gamma_x^2$ with the squared
neutral kernel. Its cost is $O(V\log V)$, not a sum over the full SD space or
an unstructured Monte Carlo estimate. More strongly, the companion
[classical comparison](classical-comparison.md) evaluates a variational trial
and its **full SD residual** in $O(V^2\log V)$ arithmetic and $O(V^2)$
memory. Its error bound is $\lVert b-Ax\rVert^2/(V\Delta)$.

Using the common rank-aware gap above, the matched $L=8,V=256$ and
$L=16,V=2048$ cases have residual error bounds approximately
$4.9026\times10^{-9}$ and $4.8801\times10^{-9}$, below the requested
$5.6981\times10^{-9}$. The observed classical times are about .0069 and .815
seconds in the final benchmark run. Those residual formulas are exact;
their displayed numerical evaluation is floating point rather than a full
interval-arithmetic end-to-end certification. The complete SD dimension is
therefore a poor measure of classical difficulty in this regime.

A classical MP2 control variate could also remove the $k=0$ term from the
quantum polynomial. Its LCU normalization becomes $L_p-1\simeq1.60478$,
instead of 2.60478, but the FFT preparation and its numerical accuracy must be
paid. That is an inherited hybrid refinement, and the deterministic classical
residual calculation is the stronger comparison for these examples. The main
ledger does not silently subtract an exactly known MP2 scalar for free.

The remaining gap to advantage is substantial: demonstrate a useful family
and precision where comparable classical scalar/residual or local expansion
methods become difficult, while preserving a certified positive gap and affordable encoding normalization and
efficiently supplied coefficients; then compile the arithmetic, phases and
fault-tolerant implementation. The present evidence supports a complete
conditional compiler and a constructive inverse refinement, not a winning
quantum application at the tested parameters.

## 8. Reproducibility and validation limits

[resource_ledger.py](../verification/resource_ledger.py) uses the Python standard library. It emits the selected
degrees, bit widths, actual small geometric CDF table, all error allocations,
confidence calculation, operation counts and control/width majorants.
[resource_results.json](../verification/results/resource_results.json) records its default volumes.
[validate_resource_ledger.py](../verification/validate_resource_ledger.py) checks shell normalization and correction
identities, exact comparator counting, fixed-point cosine evaluations against
rational Taylor/$\pi$ enclosures on selected modes, and the aggregate error
and confidence budgets at three precisions. [validation_results.json](../verification/results/resource_validation.json) records
those outcomes. These checks passed. They do not validate a complete compiled
quantum circuit, optimize its schedule, or turn the floating-point resource
report into an interval-certified hardware estimate.
