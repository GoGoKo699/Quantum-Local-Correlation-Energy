# Compiling the singles–doubles operator

[Coulomb access](coulomb-access.md) · [Gap and signal](gap-and-signal.md) · [Home](../README.md)

The SD matrix combines rank-preserving direct and exchange terms with a rank-changing block. The latter is not covered merely by a fixed-excitation-number BSE encoding.

## Normal ordering and the finite quasiparticle count

Write electron and hole creation operators as $e_x^\dagger=a_{v_x}^\dagger$ and $h_x^\dagger=a_{o_x}$. Relative to the reference, the primitive-node density fluctuation in cell $x$ is

$$
\delta n_u=
\rho_{v_x}(u)e_x^\dagger e_x
-\rho_{o_x}(u)h_x^\dagger h_x
+\rho_x(u)(e_x^\dagger h_x^\dagger+h_xe_x).
$$

The residual quartic interaction is

$$
W_N=\frac12\sum_{u,w}C_{uw}:\delta n_u\delta n_w:.
$$

Quadratic contractions are already included in the supplied Fock operator. In a sector with $m$ electrons and $m$ holes, $m\le2$, the number-density terms give at most $\binom{2m}{2}\le6$ direct pairs. In doubles these are one electron–electron pair, one hole–hole pair, and four electron–hole pairs.

The rank-preserving transition terms give

$$
\sum_{x,y}J_{xy}e_x^\dagger h_x^\dagger h_y e_y.
$$

On a distinguishable electron/hole tensor space, this is a sum of at most $m^2\le4$ lifts of

$$
\sum_{x,y}J_{xy}|x,x\rangle\langle y,y|.
$$

Each lift has norm and encoding normalization bounded by that of $J$. Subsequent antisymmetrization and valid-occupation projection do not increase the norm. Fixed-rank sorting, fermionic signs and collision masks use coordinate registers of logarithmic size.

Terms changing rank by two leave the SD subspace. The remaining number/transition cross terms form the single-to-double block $K$ and its adjoint.

## Why a target double has at most two predecessors

Consider a source single with electron at cell $r$ and hole at cell $i$. Since orbital products between distinct cells vanish, the only rank-increasing two-body route inserts a local electron–hole pair at a third cell $k$. Occupation constraints exclude $k=r$ and $k=i$.

A target double can contain at most two local electron–hole pairs. Removing one of them determines its single predecessor. Thus a target with two local pairs has at most two predecessors, one with exactly one local pair has at most one, and one with no local pair has none. Sharing only an electron or only a hole cannot introduce an additional route in this disjoint-support, one-orbital-per-band model.

For a fixed source, the allowed $k$ maps injectively to the target determinant. The retained source label and a finite route/sign label therefore suffice for reversible column generation. They are not discarded before the row/column overlap.

## Uniform source-column bound

The charge change carried by the source single is

$$
q_{r,i}=\rho_{v_r}-\rho_{o_i},
\qquad \|q_{r,i}\|_1\le2.
$$

Up to the canonical determinant sign, the target coefficient is

$$
K_{(r,i;k),(r,i)}=(\rho_k|C|q_{r,i}).
$$

With $R$ the transition-density synthesis map, the unmasked column is $R^\dagger Cq_{r,i}$. Using $CL_dC=4\pi gC$,

$$
\|R^\dagger Cq_{r,i}\|^2
\le4\pi g\gamma_*^2\langle q_{r,i},Cq_{r,i}\rangle
\le16\pi g^2\gamma_*^2G_0.
$$

Masks can only decrease this norm. If every row of a matrix has at most $r_*$ nonzeros and every column has squared norm at most $c_*^2$, Cauchy–Schwarz gives

$$
\|Kx\|^2
\le r_*\sum_j|x_j|^2\sum_i|K_{ij}|^2
\le r_*c_*^2\|x\|^2.
$$

Applying it with $r_*=2$ yields

$$
\|K\|\le4g\gamma_*\sqrt{2\pi G_0}.
$$

A charge-separated electron and hole can generate a charge–dipole field with a slower pointwise tail than neutral pair exchange. The squared column norm and bounded predecessor count, rather than a flat sparsity count, control this block.

## Constructive column and row maps

The [Green-half columns](coulomb-access.md) give a coherent column for $q_{r,i}$ with normalization $2\sqrt{gZ}$. Applying the selected block of $W^\dagger$ produces the transition response with column normalization

$$
C_{\mathrm{col}}=2\sqrt{gZ}\alpha_W.
$$

The row isometry enumerates the two candidate predecessor slots uniformly, using a disjoint failure branch when a slot is invalid. The column isometry generates the added cell and canonical target while retaining its source. Their common good labels are the ordered pair of source and target, with the determinant phase included. Use flag zero for the common good subspace, flag one for every row-map failure, and flag two for every column-map failure. Thus each failure subspace is orthogonal both to the good subspace and to the other map's failure subspace. This excludes a failure–failure overlap as well as a good–failure overlap.

Their overlap encodes

$$
K/\alpha_K,
\qquad
\alpha_K=\sqrt2C_{\mathrm{col}}
=4g\gamma_*\sqrt{2\pi Z}.
$$

The bounded predecessor list is what permits removal of the source/target bookkeeping without a volume-dependent postselection. A source norm lower bound is unnecessary.

## Combining the sectors

The off-diagonal Hermitian dilation

$$
\begin{pmatrix}0&K^\dagger\\K&0\end{pmatrix}
$$

has norm $\|K\|$ and encoding normalization $\alpha_K$. It need not be counted as two unrelated LCU terms with twice that normalization.

If $\alpha_F$ bounds the Fock encoding, a conservative full encoding is

$$
\alpha_A\le\alpha_F+6gZ+4\alpha_J+\alpha_K.
$$

All coefficients, finite route labels, local neighbor maps and fermionic signs are paid. This is a fixed-rank construction: it uses at most four quasiparticle position registers and constant route multiplicities. The count does not assert that arbitrary excitation rank can be increased without cost.

The same decomposition with the physical bound $G_0\le\overline G$ supplies the [actual gap certificate](gap-and-signal.md).
