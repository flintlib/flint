---
title: Multipoint evaluation and interpolation
summary: Evaluation of a polynomial at many points and interpolation through many points in $O(M(n) \log n)$ operations with the subproduct tree and the remainder tree, and in O(M(n)) for points in geometric progression with Bluestein's chirp transform and the method of Bostan and Schost.
complexity: $O(M(n) \log n)$ for $n$ arbitrary points (building the tree, evaluating, interpolating); $O(M(n))$ for points in geometric progression; $O(n^2)$ for the iterative methods used for short inputs.
order: 50
---

Evaluating a polynomial of degree less than $n$ at $n$ points by
Horner's rule costs $n^2$ operations, and Lagrange interpolation through
$n$ points costs the same. Both can be done in $O(M(n) \log n)$
operations by divide and conquer over a binary tree of products of
linear factors. This was found independently by Fiduccia [@Fid1972] and
by Moenck and Borodin [@MB1972; @BorMoe1974]; Strassen [@Str1973]
showed that $\Omega(n \log n)$ nonscalar multiplications are necessary
for generic points. Textbook accounts are [@GG2003, Ch. 10] and [@Ber2008],
which also describes the scaled remainder tree; Bostan, Lecerf and Schost
[@BLS2003] use the transposition principle to reduce the constants.

When the points form a geometric progression $1, q, q^2, \ldots$,
evaluation and interpolation are much cheaper: one or two
multiplications, by Bluestein's chirp transform [@Blu1970] and the
algorithm of Bostan and Schost [@BS2005].

Throughout, $x_0, \ldots, x_{n-1}$ are the points in $R = \Z/p\Z$ and
$M(n)$ is the cost of a multiplication of length $n$, with $M(n)/n$
nondecreasing.

## Generic version {#generic}

The subproduct tree is implemented once in `gr_poly`
(`_gr_poly_tree_build`, `_gr_poly_evaluate_vec_fast_precomp`,
`_gr_poly_interpolation_weights`, `_gr_poly_interpolate_fast_precomp`)
for any commutative ring. Evaluation needs only divisions by the monic
polynomials of the tree and therefore works over every ring;
interpolation needs the weights $1 / \prod_{j \ne i} (x_i - x_j)$ and
therefore points whose differences are invertible. `nmod_poly`,
`fmpz_mod_poly`, `fq_poly`, `arb_poly` and `acb_poly` call these
functions with their contexts, and multimodular interpolation over $\Z$
(`fmpz_poly_interpolate_multi_mod`) uses them modulo each prime. Over
balls, the tree algorithms give correct enclosures, but the radii can
grow much more than with evaluation point by point. The generic `gr_poly_interpolate`
currently uses Newton interpolation with $O(n^2)$ operations, which
also works over integral domains when the interpolating polynomial has
coefficients in the ring; the fast functions must be called explicitly.
The algorithms for geometric progressions are specific to `nmod_poly`.

## The subproduct tree {#tree}

Assume first $n = 2^k$. The subproduct tree has the linear polynomials
$x - x_i$ as leaves, and each internal node is the product of its two
children. Level $i$ consists of the $n/2^i$ monic polynomials

$$
P_{i,j} = \prod_{l = j 2^i}^{(j+1) 2^i - 1} (x - x_l), \qquad 0 \le j < n/2^i,
$$

of degree $2^i$, and the root is $P = P_{k,0} = \prod_l (x - x_l)$.

::: lemma {#lem-tree-cost}
The subproduct tree can be built with $O(M(n) \log n)$ operations.
:::

::: proof
Level $i + 1$ is obtained from level $i$ by $n/2^{i+1}$ products of
polynomials of degree $2^i$, costing $(n/2^{i+1}) M(2^i + 1) =
O(M(n))$ since $M(n)/n$ is nondecreasing. There are $k = \log_2 n$
levels.
:::

For general $n$ the same layout is used with $k = \lceil \log_2 n
\rceil$; the last node of each level is the product of the remaining
$n \bmod 2^i$ factors if that number is nonzero.

## Evaluation by the remainder tree {#evaluation}

::: lemma {#lem-remainder}
If $A$ divides $B$, then $(f \bmod B) \bmod A = f \bmod A$. In particular
$f(x_l) = f \bmod (x - x_l)$ can be computed by reducing $f$ successively
modulo the nodes on the path from the root to the leaf $x - x_l$.
:::

::: proof
$f - (f \bmod B)$ is a multiple of $B$, hence of $A$.
:::

::: algorithm Fast multipoint evaluation {#alg-eval}
**Input:** $f$ of length at most $n$, the subproduct tree of
$x_0, \ldots, x_{n-1}$.

**Output:** $f(x_0), \ldots, f(x_{n-1})$.

1. $r_{k,0} \gets f \bmod P$ (nothing to do if $\deg f < n$).
2. For $i = k - 1, \ldots, 0$ and all $j$: $r_{i,j} \gets r_{i+1,
   \lfloor j/2 \rfloor} \bmod P_{i,j}$.
3. Return $f(x_l) = r_{0,l}$.
:::

Correctness is [](#lem-remainder). Each level costs $n/2^i$ divisions of a
polynomial of degree less than $2^{i+1}$ by one of degree $2^i$, which
with Newton division ([](algorithms/polynomials/newton#divrem)) costs
$O(M(n))$ per level. Hence evaluation costs $O(M(n) \log n)$, in
addition to the tree. The constant can be lowered by precomputing the
inverses of the reversed tree nodes, or by the transposed and scaled
remainder trees of [@BLS2003] and [@Ber2008], which replace the
divisions by middle products; FLINT uses the plain remainder tree.

## Interpolation {#interpolation}

Let $y_0, \ldots, y_{n-1}$ be values and the $x_l$ pairwise distinct
with all $x_l - x_m$ invertible. The interpolating polynomial of degree
less than $n$ is given by Lagrange's formula

$$
f = \sum_{l=0}^{n-1} y_l w_l \frac{P}{x - x_l}, \qquad
w_l = \prod_{m \ne l} (x_l - x_m)^{-1}.
$$

::: lemma {#lem-weights}
$w_l = 1/P'(x_l)$.
:::

::: proof
$P' = \sum_l \prod_{m \ne l} (x - x_m)$, and at $x = x_l$ all terms but
the $l$-th vanish.
:::

So the weights are obtained by one multipoint evaluation of $P'$ and $n$
inversions (or one batch inversion). The sum can then be computed going
up the tree:

::: lemma {#lem-combine}
For a node $N = N_L N_R$ with children $N_L$ and $N_R$ and coefficients
$c_l$, let $S(N) = \sum_{l \in N} c_l N/(x - x_l)$, the sum over the
leaves below $N$. Then $S(N) = S(N_L) N_R + S(N_R) N_L$.
:::

::: proof
For a leaf $l$ below $N_L$, $N/(x - x_l) = N_R \cdot N_L/(x - x_l)$,
and similarly for $N_R$.
:::

::: algorithm Fast interpolation {#alg-interp}
**Input:** the subproduct tree of distinct $x_0, \ldots, x_{n-1}$, values
$y_0, \ldots, y_{n-1}$.

**Output:** the $f$ of length at most $n$ with $f(x_l) = y_l$.

1. Compute $P'$ (from the root, or from the two top nodes) and evaluate it
   at the $x_l$ by [](#alg-eval); $w_l \gets 1/P'(x_l)$.
2. $c_l \gets w_l y_l$; these are $S$ at the leaves.
3. Going up, compute $S$ at each node by [](#lem-combine).
4. Return $f = S(P)$.
:::

The cost is $O(M(n) \log n)$: one evaluation and, for the combination,
two products per node, $O(M(n))$ per level.

For few points, quadratic methods are faster:

- **Newton interpolation:** divided differences $O(n^2)$, then
  conversion from the Newton basis $\prod_{m < l} (x - x_m)$ to the
  monomial basis in $O(n^2)$.
- **Barycentric (Lagrange) interpolation:** compute $P$, the weights
  $w_l$ directly as products, and add up $y_l w_l P/(x - x_l)$, each
  quotient obtained by synthetic division (`_nmod_poly_div_root`).

## Geometric progressions {#geometric}

Let the points be $x_j = q^j$, $0 \le j < n$, and write $q = r^2$
(FLINT is given $r$, so that no square root is needed).

::: lemma Bluestein's trick {#lem-bluestein}
Let $f$ have length $\ell$. Since
$ij = \tfrac{1}{2}\bigl((i+j)^2 - i^2 - j^2\bigr)$, we have
$q^{ij} = r^{(i+j)^2} r^{-i^2} r^{-j^2}$ and

$$
f(q^j) = r^{-j^2} \sum_{i=0}^{\ell - 1} \bigl(f_i r^{-i^2}\bigr) r^{(i+j)^2}.
$$
:::

The sum is a correlation of the vector $(f_i r^{-i^2})$ with the
"chirp" sequence $c_k = r^{k^2}$, $0 \le k < \ell + n - 1$: with $a$ the
reversal of $(f_i r^{-i^2})$, $f(q^j)$ is $r^{-j^2}$ times coefficient
$\ell - 1 + j$ of $a \cdot c$, a middle product
([](algorithms/polynomials/middle-product)). So evaluation at $n$ points
of a geometric progression costs one middle product and $O(n)$
multiplications. This is the chirp z-transform [@RSR1969], also used to
compute DFTs of arbitrary length; see [@ASU1975] for its use in
polynomial evaluation.

Interpolation at a geometric progression uses a $q$-analogue of Newton
interpolation [@BS2005, §5]. With $u_i = \prod_{k=1}^{i} (q^k - 1)$ and
$Q_i = q^{i(i-1)/2}$, the Newton coefficients of the interpolant are the
low coefficients of a single product

$$
\Bigl(\sum_i \frac{y_i}{u_i} x^i\Bigr)\Bigl(\sum_i (-1)^i \frac{Q_i}{u_i} x^i\Bigr),
$$

up to scaling by $Q_i$, and the change from the Newton basis
$\prod_{m < l} (x - q^m)$ to the monomial basis is again one (transposed)
product. Hence interpolation costs two products of length $n$ and
$O(n)$ multiplications. It requires $q^k - 1$ to be invertible for
$1 \le k < n$: the points must be distinct, and for composite
modulus their differences must be invertible.

Extrapolation, computing $f(q^k)$ for $k$ outside $[0, \ell)$ from the
values $f(q^0), \ldots, f(q^{\ell - 1})$ of a polynomial of length $\ell$,
is a single middle product as well:

::: lemma {#lem-extrapolate}
For $k \notin [0, \ell)$,

$$
f(q^k) = \prod_{j=0}^{\ell-1} (q^k - q^j) \sum_{i=0}^{\ell-1}
\frac{f(q^i) \, w_i \, q^{-i}}{q^{k-i} - 1}, \qquad
w_i = \prod_{j \ne i} (q^i - q^j)^{-1}.
$$

The sum is a convolution of the sequence $f(q^i) w_i q^{-i}$ with the
sequence $1/(q^m - 1)$, and the prefactors are products of the numbers
$q^m - 1$ times powers of $q$.
:::

::: proof
This is Lagrange's formula $f(z) = P(z) \sum_i f(q^i) w_i/(z - q^i)$,
valid for $z$ not a node, with $z = q^k$ and $q^k - q^i = q^i (q^{k-i} - 1)$.
:::

## Use in FLINT {#use}

**Arbitrary points.**

- `_nmod_poly_tree_alloc`, `_nmod_poly_tree_build` and
  `_nmod_poly_tree_free` manage the subproduct tree. `tree[i]` holds the
  nodes of level $i$ consecutively (each of length $2^i + 1$, the last
  possibly shorter). Level 1 is computed directly from the points and
  the root is not built, since evaluation reduces modulo the highest
  level below the length of $f$ and interpolation forms the top product
  only once, to obtain $P'$.
- `nmod_poly_evaluate_nmod_vec` uses Horner's rule at each point for
  polynomials of length less than 32 (`nmod_poly_evaluate_nmod_vec_iter`)
  and [](#alg-eval) otherwise (`nmod_poly_evaluate_nmod_vec_fast`), which
  builds the tree for each call. `_nmod_poly_evaluate_nmod_vec_fast_precomp`
  takes a prebuilt tree, and accepts polynomials longer or shorter than
  the number of points: it first reduces $f$ modulo all nodes of the
  highest built level whose nodes have degree at most $\deg f$. The divisions use
  `_gr_poly_divrem` with the `nmod` context, so they switch to Newton
  division at the cutoffs of [](algorithms/polynomials/newton#divrem).
- `nmod_poly_interpolate_nmod_vec` uses Newton interpolation for $n < 6$
  (`nmod_poly_interpolate_nmod_vec_newton`), the barycentric form for
  $n < 16$ (`nmod_poly_interpolate_nmod_vec_barycentric`) and
  [](#alg-interp) otherwise (`nmod_poly_interpolate_nmod_vec_fast`).
  `_nmod_poly_interpolation_weights` and
  `_nmod_poly_interpolate_nmod_vec_fast_precomp` split the fast method
  into its two steps, so that the weights can be reused for many sets of
  values at the same points.
- `nmod_poly_product_roots_nmod_vec` computes $P$ alone, by multiplying
  in the linear factors one at a time for $n < 20$ and by splitting the
  points into two halves recursively otherwise.

**Geometric progressions.** `nmod_geometric_progression_init` prepares
the evaluation, interpolation and extrapolation data for the points
$1, r^2, r^4, \ldots, r^{2(\ell - 1)}$; `_nmod_geometric_progression_init_function`
selects which of the three to prepare and can also store `fft_small`
transforms of the fixed factors (the chirp sequence, the interpolation
series and the extrapolation series), so that each later product needs
two transforms instead of three. The products of evaluation and
extrapolation are middle products in which the fixed factor $v$ is the
long operand, and these are computed as transposed products
([](algorithms/polynomials/middle-product#transposed-fft)): the
precomputation stores the images of $v$ through the transposed inverse
transform (`_mul_precomp_init_tr`), one transform of length $v_n$
rounded up to a multiple of 256 per prime, with $v_n = 2\ell - 1$ for
the chirp sequence and $v_n = \ell - 1$ for each of the two extrapolation
series, instead of a power of two at least $v_n$; each call then
transforms the reversed short operand and applies the transposed forward
transform (`_nmod_poly_mul_tr_precomp`). The interpolation series is
stored as an ordinary transform (`_mul_precomp_init`). A call uses the
stored images of $v$ only if the product would use `fft_small` anyway,
if it fits, and if this is estimated to be cheaper than what
`_nmod_poly_mulmid_fft_small` would do without the precomputation. For
$p < 2^{\bits - 1}$ the
evaluation and interpolation data are computed with [Shoup multiplication](algorithms/arithmetic/shoup)
for the many products by powers of $r$. Then

- `nmod_poly_evaluate_geometric_nmod_vec_fast` evaluates by
  [](#lem-bluestein) with one middle product
  (`nmod_poly_evaluate_geometric_nmod_vec_iter` uses Horner's rule);
- `nmod_poly_interpolate_geometric_nmod_vec_fast` interpolates by the
  method of Bostan and Schost;
- `nmod_poly_extrapolate_geometric` extrapolates by [](#lem-extrapolate)
  with one middle product.

These are the methods of choice whenever the evaluation points can be
chosen freely, for example in evaluation–interpolation schemes for
products of multivariate or matrix polynomials.

The functions are described, with examples, in [](nmod_poly/multipoint).

::: functions _nmod_poly_tree_build _nmod_poly_evaluate_nmod_vec_fast_precomp _nmod_poly_interpolate_nmod_vec_fast_precomp nmod_poly_evaluate_nmod_vec nmod_poly_interpolate_nmod_vec nmod_poly_evaluate_geometric_nmod_vec_fast nmod_poly_interpolate_geometric_nmod_vec_fast nmod_poly_extrapolate_geometric :::
