---
title: Multiplying matrices over the integers
kind: topic
summary: The methods FLINT uses to multiply integer matrices — word-size dot products, Waksman and Strassen–Winograd on multiprecision entries, multimodular products through nmod_mat, floating-point gemm modulo small primes, and transforming each entry once with the small-prime FFT — their costs, and when fmpz_mat_mul uses which.
order: 35
---

Multiplying an $m \times k$ by a $k \times n$ integer matrix with entries
of $a$ and $b$ bits produces entries of up to
$c = a + b + \log_2 k$ bits. Two parameters decide the best method: the
dimension $d = \min(m, k, n)$ and the entry size. For small entries the
cost is dominated by memory traffic and by the overhead of the `fmpz`
representation; for large matrices of moderate entries, by the $mkn$
multiply-adds; for large entries, by the cost of each integer
multiplication. `fmpz_mat_mul` implements the choice; the exact tests
are in [](fmpz_mat/mul#dispatch).

The generic `gr_mat_mul` over the integers (`gr_ctx_init_fmpz`) calls
`fmpz_mat_mul`. Conversely, several other modules reduce their matrix
products to it: `fmpq_mat_mul` after clearing denominators of rows and
columns, `fmpz_mod_mat_mul` before reducing, `fmpz_poly_mat_mul_KS` and
`nmod_poly_mat_mul_KS` after Kronecker substitution, `arb_mat_mul_block`
and `nfloat`'s matrix multiplication on blocks of exponent-aligned
mantissas, and `padic_mat_mul`.

## The methods {#methods}

| Method | Function | Cost | Best for |
|---|---|---|---|
| Dot products in machine words | `_fmpz_mat_mul_small_internal`, `_fmpz_mat_mul_double_word_internal` | $mkn$ word multiply-adds, 1 to 5 words | entries up to two words, $d$ up to a few hundred |
| Classical dot products | `fmpz_mat_mul_classical` | $mkn$ multiplications of entries | tiny $k$, irregular entries |
| Waksman | `fmpz_mat_mul_waksman` | $\approx mkn/2$ multiplications, $O(mkn)$ additions | $d < 20$, balanced entries of thousands of bits |
| Strassen–Winograd | `fmpz_mat_mul_strassen` | $O(d^{\log_2 7})$ multiplications of entries | $d \ge 8$, entries of hundreds of bits |
| Floating-point `gemm` | `fmpz_mat_mul_blas` | $\approx (c / 20)$ `dgemm` calls, or one if $c \le 53$ | large $d$, $c$ up to a few hundred bits |
| Multimodular | `fmpz_mat_mul_multi_mod` | $\approx c/59$ products over $\Z/p\Z$, plus reductions | large $d$, entries of up to thousands of bits |
| Transformed entries | `fmpz_mat_mul_fft_small` | $mk + kn + mn$ transforms, $mkn$ pointwise products | entries of many thousands of bits |
| Schönhage–Strassen entries | `fmpz_mat_mul_fft` | as above, with the `fft` module | not used by `fmpz_mat_mul` |

### Machine words {#words}

When all entries fit in a signed word minus two bits (small `fmpz`), the
product can be computed without touching the multiprecision machinery:
the kernels read the words, multiply them with one hardware
multiplication, and accumulate in one, two or three words depending on
the bound $c$. `_fmpz_mat_mul_small_internal` transposes $B$ into an
array of words and processes blocks of 16 rows of $A$ against 128 entries
of the inner dimension, which keeps the working set in cache, and splits
the rows over threads. Entries of up to two words are handled likewise by
`_fmpz_mat_mul_double_word_internal`, which multiplies two-word integers
with four word products and accumulates in four or five words. For
dimensions above about 1000 these lose to Strassen–Winograd (when $c$ fits
in a word) or to multimodular multiplication.

### Fewer multiplications {#fewer}

For entries of thousands of bits, each multiplication of entries costs
$M(b)$, much more than an addition, so algorithms that trade
multiplications for additions win for small dimensions. Waksman's
algorithm [@Waksman1970] computes, for each pair of inner indices, the
products $(a_{i,2l} + b_{2l+1,j})(a_{i,2l+1} + b_{2l,j})$, whose sum over
$l$ is the dot product plus terms depending only on row $i$ or column
$j$; it needs $mkn/2 + O(mk + kn)$ multiplications. Strassen–Winograd
[@Str1969; @Win1971], in the operation order of Bodrato
[@Bodrato2010], needs 7 instead of 8 half-size products per level at the
cost of 15 additions; FLINT applies one level and recurses through
`fmpz_mat_mul`. See [](algorithms/linear-algebra/strassen-winograd).

### Products modulo primes {#multimodular}

The entries of $C$ are less than $2^c$ in absolute value, so $C$ is
determined by $C$ modulo moduli $p_1, \ldots, p_r$ with product at least
$2^{c + 1}$ (or $2^c$ if all entries are nonnegative). The multimodular
method reduces $A$ and $B$ modulo each $p_i$, multiplies with
`nmod_mat_mul`, and reconstructs $C$ by the Chinese remainder theorem
([](algorithms/linear-algebra/multimodular#crt)). FLINT uses the modulus
$2^{63}$, for which reduction is free, and primes above $2^{59}$, so
$r \approx 1 + (c - 63)/59$; if $c \le 63$, the single modulus $2^c$.

::: proposition {#prop-multimod-cost}
For $d \times d$ matrices with entries of $b$ bits, multimodular
multiplication costs
$$
O\big(r \cdot d^\omega + d^2 \cdot T(b)\big), \qquad r = O(b/\beta + \log d / \beta + 1),
$$
word operations, where $d^\omega$ is the cost of `nmod_mat_mul` and
$T(b) = O(M(b) \log(b/\beta))$ that of reducing a $b$-bit integer modulo
$r$ primes, or of reconstructing it, with a subproduct tree.
:::

::: proof
There are $2d^2$ entries to reduce and $d^2$ to reconstruct, each by a
tree of products of the $r$ moduli [@GG2003, §10.3], and $r$ matrix
products modulo word-size moduli.
:::

Compared with the $d^3 M(b)$ of the classical algorithm, the matrix
products cost $d^\omega$ word operations per prime, that is, about
$d^\omega b/\beta$ in total, which for $b$ up to a few thousand bits is
much less than $d^3 M(b)$, and the conversions are only quadratic in $d$.
Since `nmod_mat_mul` runs at nearly the floating-point peak of the
processor ([](algorithms/linear-algebra/matrix-multiplication)), this is
the method of choice for large matrices.

### Floating point {#floating-point}

`fmpz_mat_mul_blas` is a multimodular method with primes small enough for
a `double` matrix product to be exact: residues in the symmetric range
below $p/2$ with $k (p/2)^2 < 2^{53}$, so that every dot product is
computed exactly ([](algorithms/linear-algebra/gemm#exactness)). That
gives primes of about $26 - \frac12 \log_2 k$ bits (at most about 21 bits
when $c > 200$, to keep the reconstruction cheap), and about
$c/(26 - \frac12 \log_2 k)$ calls to `flint_dgemm`. If
$a + b + \log_2 k \le 53$, a single `dgemm` on the entries suffices. With
a tuned BLAS library or FLINT's own kernels
([](algorithms/linear-algebra/gemm)) this is the fastest method for large
matrices with entries up to a few hundred bits: `fmpz_mat_mul` uses it
for $d > 50$ when $c \le 53$, and from about $d = 250$ to $1300$
depending on $c$.

### Transformed entries {#transformed}

For entries of $b \gg 1000$ bits, each product of entries costs a
multiplication of large integers, which with an FFT means two forward
transforms, a pointwise product and an inverse transform. In a matrix
product each entry of $A$ takes part in $n$ products and each entry of
$B$ in $m$, so transforming every entry once and adding the $k$
pointwise products of each entry of $C$ in the transformed domain saves
almost all transforms.

::: proposition {#prop-transformed}
With transforms of length $L = O(b/\log b)$ points, multiplying through
transformed entries costs
$$
O\big((mk + kn + mn) \, L \log L + mkn \, L\big)
$$
operations, against $O(mkn \, L \log L)$ for multiplying the entries
separately.
:::

::: proof
Each of the $mk + kn$ input entries and $mn$ output entries is
transformed once, and each of the $mkn$ products is a pointwise product
of length $L$. The transform length must accommodate the accumulated sum
of $k$ products, which adds $\log_2 k$ bits to the coefficients and
changes $L$ by a constant factor.
:::

`fmpz_mat_mul_fft_small` implements this with the transformed integer
ring of `gr`, built on the small-prime FFT of [](fft_small)
([](algorithms/fft/multimodular#reuse)). The transforms of the inputs
must be stored, $(mk + kn) L$ numbers, which the function limits with
`flint_fft_small_max_transformed_ring_size`, falling back to row and
column blocks and, if necessary, to blocks of the inner dimension.
`fmpz_mat_mul` uses it when the smaller entry size exceeds about
$7000 + 0.15 d^2$ bits (10000 for $d \le 5$). The older
`fmpz_mat_mul_fft` does the same with the Schönhage–Strassen transform of
the [](fft) module, for entries of at least 2000 bits.

## Squaring and powers {#squaring}

For $A^2$, the transformed method keeps one pool of transforms, and
Bodrato's variant of Strassen–Winograd [@Bodrato2010] uses squarings of
half-size blocks where possible. `fmpz_mat_sqr` applies Bodrato's
formulas for dimensions up to 3 (and 4 with large entries) and otherwise
`fmpz_mat_mul`. `fmpz_mat_pow` uses binary powering.

## Rational and polynomial matrices {#related}

- `fmpq_mat_mul` multiplies $\operatorname{diag}(r)^{-1} A'$ and
  $B' \operatorname{diag}(s)^{-1}$, where $A'$ and $B'$ are $A$ and $B$
  with rows, respectively columns, multiplied by their least common
  denominators $r_i$ and $s_j$, so that one integer product suffices.
- Polynomial matrices over $\Z$ and $\Z/n\Z$ are multiplied by Kronecker
  substitution into integer matrices, or by evaluation and interpolation;
  see [](algorithms/linear-algebra/polynomial-matrices).
