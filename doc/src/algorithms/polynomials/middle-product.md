---
title: Middle products and the transposition principle
summary: The middle product, a window of coefficients of a product, is the transpose of polynomial multiplication and costs no more than an ordinary product of the corresponding size; FLINT computes it with dot products, a cyclic convolution whose wrap-around misses the window, or transposed FFTs, and uses it in Newton iteration and geometric evaluation.
complexity: A middle product of lengths $m + n - 1$ and $n$ (output length $m$) costs $mn$ multiplications classically and about one product of lengths $m$ and $n$ with the FFT.
order: 30
---

Many algorithms need only a window of the coefficients of a product. The
most important case is the *middle product*: given $a$ of length
$m + n - 1$ and $b$ of length $n$, the $m$ coefficients of $ab$ of
index $n - 1$ to $m + n - 2$. These are exactly the coefficients in
which every coefficient of $b$ meets a coefficient of $a$. The middle
product is the transpose of multiplication, and the transposition
principle says that it can be computed as fast as a product of lengths
$m$ and $n$, although the full product $ab$ has length $m + 2n - 2$.
The middle product was named and systematically studied by Hanrot,
Quercia and Zimmermann [@HQZ2004]. The transposition principle comes
from Tellegen's work on electrical networks and was brought into
algebraic complexity by Fiduccia [@BCS1997, Ch. 13]; Shoup [@Sho1995]
and Bostan, Lecerf and Schost [@BLS2003] made systematic use of it for
polynomial algorithms.

FLINT generalises the notion: `_nmod_poly_mulmid` computes coefficients
$m_1 \le k < m_2$ of a product of any lengths, so that `mullow` is the
case $m_1 = 0$ and `mul` the case $m_1 = 0$, $m_2 = \ell_1 + \ell_2 - 1$.

## Generic version {#generic}

`gr_poly_mulmid` computes the window $[m_1, m_2)$ of a product over any
ring through the ring method `GR_METHOD_POLY_MULMID`. Rings with a fast
multiplication supply a middle product with it: `nmod` uses
`_nmod_poly_mulmid` as described below, and `fmpz`, `fmpq`, `fmpz_mod`,
`mpn_mod`, the `fq` family, the ball rings, `fmpzi` and polynomial rings
(by bivariate Kronecker substitution) have their own methods. The
default `_gr_poly_mulmid_generic` computes the low product of length
$m_2$ with the ring's `mullow` method and discards the first $m_1$
coefficients when that method is not the classical default, and
otherwise uses the classical middle product, which computes only the
requested coefficients with one dot product each. The generic Newton
iterations ([](algorithms/polynomials/newton#generic)) check whether the
ring has its own middle product and use it only in that case, since
the fallback through a low product gains nothing.

## The transposition principle {#tellegen}

A *linear straight-line program* computes $y = Ax$ for a fixed matrix
$A \in K^{r \times s}$ by a sequence of additions, subtractions and
multiplications by constants applied to the inputs $x_1, \ldots, x_s$
and to earlier results. Its length is the number of such operations.

::: theorem Tellegen's principle {#thm-tellegen}
If $A \in K^{r \times s}$ has no zero rows and no zero columns and $y = Ax$
is computed by a linear straight-line program of length $L$, then
$x' = A^T y'$ is computed by a linear straight-line program of length
$L + r - s$.
:::

::: proof
(Sketch; see [@BCS1997, Ch. 13] for details.) Represent the program as
a directed acyclic graph whose vertices are the inputs, the
intermediate results and the outputs, with an edge, weighted by the
constant factor, from each operand of an operation to its result. Then
$A_{ij}$ is the sum, over all paths from input $j$ to output $i$, of the
products of the weights. Reversing every edge gives a graph in which the
same sums describe $A^T$; the multiplications by constants stay on
their edges. A vertex with $d$ incoming edges costs $d - 1$ additions,
so with $E$ edges and $V$ vertices the original program has
$\sum (\mathrm{indeg}(v) - 1) = E - V + s$ additions, the sum running
over all vertices except the $s$ inputs. In the reversed graph the
in-degrees are the original out-degrees and the sources are the $r$
outputs, so it has $E - V + r$ additions. The hypothesis on $A$ ensures
that no input or output is isolated.
:::

The principle turns any algorithm for a linear map into one for its
transpose with essentially the same cost. Transposing an algorithm by
hand is mechanical but error-prone; [@BLS2003] give a methodology and
many examples.

## Middle product as transposed multiplication {#transpose}

For a fixed $b$ of length $n$, multiplication by $b$ is the linear map
$\mu_b: K^m \to K^{m+n-1}$, $a \mapsto ab$, whose matrix
$(b_{k-i})_{k,i}$ is a Toeplitz matrix.

::: lemma {#lem-transpose}
$\mu_b^T(c) = \bigl( \sum_{j=0}^{n-1} c_{i+j} b_j \bigr)_{0 \le i < m}$
for $c$ of length $m + n - 1$, and these are the coefficients
$n - 1, \ldots, m + n - 2$ of $c \cdot \operatorname{rev}(b)$, where
$\operatorname{rev}(b) = x^{n-1} b(1/x)$.
:::

::: proof
The entry of $\mu_b$ in row $k$ and column $i$ is $b_{k-i}$ (zero
outside $0 \le k - i < n$), so $(\mu_b^T c)_i = \sum_k b_{k-i} c_k =
\sum_j b_j c_{i+j}$. On the other hand,
$(c \operatorname{rev}(b))_{i+n-1} = \sum_l c_{i+n-1-l} b_{n-1-l}$,
which is the same sum with $j = n - 1 - l$.
:::

By [](#thm-tellegen), any algorithm that multiplies polynomials of
lengths $m$ and $n$ yields an algorithm for the middle product with
essentially the same number of operations. [@HQZ2004] give explicit
transposed versions of Karatsuba's and the FFT algorithms; for the FFT
the transposed algorithm is particularly simple, as shown next.

## Classical and FFT middle products {#fft}

The classical middle product computes each requested coefficient as a
dot product, $mn$ multiplications in total for the HQZ middle product.
In FLINT these dot products use the delayed reduction of
[](algorithms/linear-algebra/dot-products) (`_nmod_poly_mulmid_classical`).

An FFT computes products modulo $x^N - 1$, that is cyclic convolutions,
and a middle product only needs a convolution length large enough for
the wrap-around to miss the requested window.

::: lemma Wrap-around {#lem-wrap}
Let $c = ab$ have length $z$ and let $0 \le m_1 < m_2 \le z$. If
$m_2 \le N$ and $z - N \le m_1$, then coefficients $m_1$ to $m_2 - 1$ of
$ab \bmod (x^N - 1)$ equal those of $c$.
:::

::: proof
Coefficient $k < N$ of $c \bmod (x^N - 1)$ is $\sum_{t \ge 0} c_{k + tN}$.
For $m_1 \le k < m_2$ we have $k < N$ and $k + N \ge m_1 + N \ge z$, so
only the term $t = 0$ is nonzero.
:::

For the HQZ middle product, $z = m + 2n - 2$, $m_1 = n - 1$ and
$m_2 = m + n - 1$, so any $N$ with $m + n - 1 \le N \le m + 2n - 2$
works: the transforms then have the length of the *longer operand*, as
for a product of lengths $m$ and $n$, instead of the length $m + 2n - 2$
of the full product. A radix-2 FFT needs $N$ to be a power of two, which
need not exist in that range (for $m = n$ the range is
$[2n - 1, 3n - 2]$); a truncated transform of length below $2^d$ does not
compute a product modulo $x^N - 1$. The transposed FFT
([](#transposed-fft)) has no such restriction.

::: algorithm FFT middle product {#alg-fft-mulmid}
**Input:** $a$ of length $\ell_1$, $b$ of length $\ell_2$, a window
$[m_1, m_2)$ of the product, of length $z = \ell_1 + \ell_2 - 1$.

**Output:** coefficients $m_1, \ldots, m_2 - 1$ of $ab$.

1. If there is a power of two $N$ with $\max(\ell_1, \ell_2, m_2) \le N
   \le z \le m_1 + N$, take it. Otherwise take a transform length
   $N \ge z$ (no wrap-around).
2. Compute the forward transforms of $a$ and $b$ of length $N$, multiply
   pointwise, and transform back.
3. Return coefficients $m_1$ to $m_2 - 1$ of the result.
:::

Correctness is [](#lem-wrap). This is what `fft_small` does in
`_nmod_poly_mul_mid_mpn_ctx` (via `_fft_small_plan_set_window_nmod`),
with truncated transforms rounded to blocks of 256 when there is no
suitable power of two; [](algorithms/fft/multimodular#windows) treats
the same wrap-around for windows of integer products. Each transform is computed modulo one to four
primes and the results are combined by Chinese remaindering, as for
ordinary products ([](algorithms/polynomials/multiplication#fft-small)).

### Transposed FFT {#transposed-fft}

The transposition principle can also be applied to the FFT product
itself. Fix $b$ of length $n$ and an output length $m$, and let
$T \ge m + n - 1$ be the truncation of a transform of length $2^d$.
Multiplication by $b$ is computed as
$$\mu_b = R \circ I \circ D \circ F,$$
where $F : K^m \to K^T$ evaluates a polynomial of length $m$ at the
first $T$ points (a truncated forward transform), $D$ multiplies
pointwise by the transform of $b$ and a normalising constant, $I$ is the
truncated inverse transform of length $T$, and $R : K^T \to K^{m+n-1}$
keeps the first $m + n - 1$ coefficients.

::: theorem Transposed convolution {#thm-transposed-fft}
With this notation, $\mu_b^T = F^T \circ D \circ I^T \circ R^T$, where
$R^T$ pads a vector of length $m + n - 1$ with zeros to length $T$. Hence
$v \mapsto F^T(D(I^T(v)))$ computes the middle product
$\bigl(\sum_j b_j v_{i+j}\bigr)_{0 \le i < m}$ of [](#lem-transpose)
with the transforms of a product of lengths $m$ and $n$, and the
transposed transforms $F^T$ and $I^T$ cost as much as $F$ and $I$.
:::

::: proof
The factorisation of $\mu_b$ holds as an identity of linear maps,
because the truncated transforms compute the product exactly for every
input of length $m$ once $T \ge m + n - 1$. Transposing reverses the
order of the factors, and $D$ is diagonal. The transforms are linear
straight-line programs (butterflies, copies and zeroings), so by
[](#thm-tellegen) their transposes are computed by programs of the same
length; in fact each butterfly is replaced by its transpose and the
order of the butterflies is reversed, which gives the same numbers of
additions and multiplications.
:::

::: algorithm Transposed FFT middle product {#alg-tr-fft}
**Input:** $b$ of length $n$, $v$ of length at most $m + n - 1$.

**Output:** $z_i = \sum_{j < n} b_j v_{i+j}$ for $0 \le i < m$ (with $v$
padded with zeros); for $b = \operatorname{rev}(f)$ these are
coefficients $n - 1$ to $m + n - 2$ of $fv$.

1. Choose the transform length $2^d$ and truncation $T$ of a product of
   lengths $m$ and $n$.
2. Compute the forward transform $\hat b$ of $b$ (input truncation $n$,
   output truncation $T$).
3. Pad $v$ to length $T$ and apply the transposed inverse transform.
4. Multiply pointwise by $\hat b$ and the normalising constant.
5. Apply the transposed forward transform with input truncation $m$ (it
   reads $T$ values and writes $m$) and return its first $m$ entries.
:::

Correctness is [](#thm-transposed-fft). Compared with
[](#alg-fft-mulmid), the transform length is that of a product of
lengths $m$ and $n$ for every $m$ and $n$, not only when a power of two
lies in $[m + n - 1, m + 2n - 2]$. No wrap-around is used: the transpose
of a cyclic convolution is a cyclic correlation, which is not the
middle product. The same idea gives transposed versions of every
FFT-based multiplication, see [@BLS2003] and, for the truncated Fourier
transform, [@vdH2004].

In FLINT, `sd_fft_trunc_tr` and `sd_ifft_trunc_tr` are the transposes of
the truncated transforms of `fft_small`
([](algorithms/fft/small-prime-fft#transposed)), and
`_nmod_poly_mul_tr_mpn_ctx` is [](#alg-tr-fft) modulo one to four primes
with Chinese remaindering. `_nmod_poly_mulmid_fft_small` computes a
window as a window of the full product with [](#alg-fft-mulmid) or as a
transposed product, whichever has the strictly shorter transform; the
former wins only when its wrap-around reaches the same length, where it
is a few percent faster. For the balanced middle product ($m = n$, as in
Newton iteration) without a power of two in $[2n - 1, 3n - 2]$, the
transform shrinks from about $3n$ to $2n$ points.

The same reduction applies to any window $[m_1, m_2)$ of a product of
lengths $\ell_1$ and $\ell_2$. First the coefficients of the operands
that no output coefficient uses are dropped: those of index at least
$m_2$, and the low coefficients of each operand that meet the other one
only below $m_1$. Let then $s$ of length $n$ be the shorter operand, $l$
of length $\ell$ the longer one and $m = m_2 - m_1$. Coefficient
$m_1 + i$ of $sl$ is $\sum_j \operatorname{rev}(s)_j v_{i+j}$ with
$v_t = l_{t + m_1 - (n - 1)}$ (zero outside $[0, \ell)$), a transposed
product of lengths $n$ and $m$. When $m < \ell$ this costs a product of
the shorter operand and the window instead of a product of both
operands; since the operands enter symmetrically, they may be given in
either order. `_nmod_poly_mulmid_fft_small`, the integer windows of
`_mpn_ctx_mpn_mul_range` (on the polynomials of chunks, see
[](fft_small/mpn_mul#algorithm)) and `radix_mulmid_fft_small` use it
whenever $m < \ell$ and its transform is strictly shorter than that of
[](#alg-fft-mulmid). The integer versions skip it for squarings, and the
`nmod_poly` version also requires that the window does not extend beyond
the full sums on both sides.

The transposed product $z = F^T(\hat b \cdot I^T(v))$ is bilinear in $b$
and $v$, so either factor can be prepared in advance. When $v$ is fixed
and many short $b$ are applied to it, storing $I^T(v)$ leaves two
transforms per product, $F(b)$ and $F^T$, against three for the plain
transposed product, and also for [](#alg-fft-mulmid) with a stored
transform of $v$. The transform length then depends only on
$v_n \ge m + n - 1$ (in `fft_small`, $v_n$ rounded up to a multiple of
256, instead of a power of two or the length $v_n + n - 1$ of the full
product). This is `_mul_precomp_init_tr` with `_nmod_poly_mul_tr_precomp`,
used for the middle products of evaluation and extrapolation at
geometric progressions
([](algorithms/polynomials/subproduct-tree#geometric)).

Kronecker substitution has no such shortcut in FLINT:
`_nmod_poly_mulmid_KS` drops the operand coefficients that cannot
contribute to the window and then forms the whole integer product.

## Uses {#uses}

**Newton iteration.** When $g$ approximates $1/Q$ to $m$ terms, the
first $m$ coefficients of $Qg$ are known in advance ($1, 0, \ldots, 0$),
and the next step needs only coefficients $m$ to $2m - 1$. With $Q$
truncated to length $2m$ this is a middle product with output length
$m$, costing about one product of length $m$ instead of a product of
lengths $2m$ and $m$. The same saving appears in division, square roots
and exponentials; see [](algorithms/polynomials/newton).

**Geometric progressions.** Evaluating a polynomial at the points
$1, q, q^2, \ldots$ reduces, by Bluestein's chirp trick, to a single
middle product ([](algorithms/polynomials/subproduct-tree#geometric)).

**Transposed algorithms.** Transposing multipoint evaluation gives
the "transposed remainder tree" of [@BLS2003] and Bernstein's scaled
remainder tree [@Ber2008], and transposing modular multiplication gives
Shoup's power projection [@Sho1995], used in minimal polynomial and
factoring algorithms. FLINT's `nmod_poly` does not currently use these
transposed variants: multipoint evaluation uses the ordinary remainder
tree ([](algorithms/polynomials/subproduct-tree)).

## Use in FLINT {#use}

- `_nmod_poly_mulmid` and `nmod_poly_mulmid` compute any window
  $[m_1, m_2)$ of a product, choosing between the classical method,
  Kronecker substitution, `_nmod_poly_mul` (for the full window) and the
  FFT, see [](algorithms/polynomials/multiplication#cutoffs). The forced
  variants are `_nmod_poly_mulmid_classical`, `_nmod_poly_mulmid_KS` and
  `_nmod_poly_mulmid_fft_small`; the last one chooses between
  [](#alg-fft-mulmid) (`_nmod_poly_mul_mid_default_mpn_ctx`) and
  [](#alg-tr-fft) (`_nmod_poly_mul_tr_default_mpn_ctx`).
- In `fft_small`, `fft_small_ifft_tr_nmod` and `fft_small_fft_tr` give
  the transposed transforms on transformed operands
  ([](fft_small/transform-space#transposed)), and `_mul_precomp_init_tr`
  with `_nmod_poly_mul_tr_precomp` the transposed product with a stored
  $I^T(v)$. The integer middle products `flint_mpn_mulmid` (through
  `_mpn_ctx_mpn_mul_range`) and `radix_mulmid` (through
  `radix_mulmid_fft_small`) use [](#alg-tr-fft) for short windows as
  described above; the `fmpz_poly` and `mpn_mod` middle products of
  `fft_small` use only [](#alg-fft-mulmid).
- The `gr` context for $\Z/n\Z$ installs `_nmod_poly_mulmid` as its
  middle product method. The generic Newton iterations of `gr_poly`
  (`_gr_poly_inv_series_newton`, `_gr_poly_div_series_newton`,
  `_gr_poly_rsqrt_series_newton`, `_gr_poly_sqrt_series_newton`,
  `_gr_poly_exp_series_newton`), which `nmod_poly` calls for long series,
  use a middle product whenever the ring provides one, and a low
  product followed by discarding coefficients otherwise.
- `nmod_poly_evaluate_geometric_nmod_vec_fast` and
  `nmod_poly_extrapolate_geometric` are each one middle product plus
  linear work; with a precomputation made with bit 3 of `function` set
  (`_nmod_geometric_progression_init_function`), these middle products
  use a stored $I^T(v)$.
- `fmpz_poly_mulmid` provides the same operation over $\Z$.

::: functions _nmod_poly_mulmid nmod_poly_mulmid _nmod_poly_mulmid_classical _nmod_poly_mulmid_KS _nmod_poly_mulmid_fft_small _nmod_poly_mul_tr_mpn_ctx :::
