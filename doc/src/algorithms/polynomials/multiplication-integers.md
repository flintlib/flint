---
title: Multiplying polynomials over Z
kind: topic
summary: How fmpz_poly multiplies, depending on the length and the coefficient size - word-size accumulation, classical and packed-limb products, Karatsuba with Toom-Cook interpolation in the polynomial variable, Kronecker substitution, the multimodular small-prime FFT, Schönhage-Strassen convolutions and classical multiplication in transform space - with the exactness lemmas and the dispatcher that chooses among them.
order: 5
---

A polynomial over $\Z$ has two sizes: its length $n$ and the bit size
$b$ of its coefficients. A product of two polynomials of length $n$
with $b$-bit coefficients has $2n - 1$ coefficients of at most $2b +
\log_2 n$ bits, and its cost is bounded below by the cost of
multiplying integers of about $nb$ bits. Every algorithm described below
reaches that cost in some range and is far from it in others: for $n =
3$ and $b = 10^6$ the problem is a handful of huge integer products,
for $n = 10^6$ and $b = 10$ it is essentially a product over a small
finite field, and in between both dimensions matter. FLINT therefore
selects among more algorithms for $\Z[x]$ than for $\Z/n\Z[x]$
([](algorithms/polynomials/multiplication)), in a single dispatcher,
`_fmpz_poly_mulmid`, which computes an arbitrary window $[n_{lo},
n_{hi})$ of coefficients of the product and serves full, low, high and
middle products and squares. The generic `gr_poly` functions over the
ring $\Z$ (`gr_ctx_init_fmpz`) call the same dispatcher through the
method table of the ring, so generic algorithms such as power series
composition or reversion over $\Z$ inherit it.

Throughout, $a$ and $b$ have lengths $\ell_1 \ge \ell_2$ and
coefficients of at most $\beta_1$ and $\beta_2$ bits, $\beta = \beta_1
+ \beta_2$, and $c = ab$. Each coefficient of $c$ is a sum of at most
$\ell_2$ products, so

$$
|c_k| \le \ell_2\, 2^{\beta} \qquad (0 \le k < \ell_1 + \ell_2 - 1).
$$

## Direct algorithms {#direct}

**Word accumulation.** If $\beta_1, \beta_2 \le 62$ (small `fmpz`
values) and $\beta + \lceil \log_2(\ell_2 + 1) \rceil$ fits in one
signed word, the products are accumulated directly in machine words,
with loops that the compiler vectorises; if it fits in two words, each
coefficient is accumulated in a double word with `smul_ppmm` and
`add_ssaaaa`. No `fmpz` arithmetic happens until the result is
written.

**Classical multiplication.** Each output coefficient is computed as
one dot product of `fmpz` values with `_fmpz_vec_dot_general`;
squares use $c_k = 2 \sum_{i < k/2} a_i a_{k-i} +
a_{k/2}^2$. The cost is the number of coefficient products in the
window, which makes it the best method for very short operands and for
narrow windows.

**Packed limbs.** `_fmpz_poly_mulmid_mpn` copies the coefficients into
fixed-size limb arrays (unsigned, two's complement, sign-magnitude or
with a bias, chosen by `_flint_mpn_poly_mulmid_params`) and computes
each output coefficient modulo $2^{\bits s}$ for an $s$ large enough
for the exact value, with classical, Karatsuba and Karatsuba
middle-product algorithms on the packed representation. This removes
the overhead of the `fmpz` type, which dominates for coefficients of a
few limbs.

## Karatsuba and Toom–Cook in the polynomial variable {#toom}

For large coefficients, the cost is the number of coefficient products,
and evaluation–interpolation in the polynomial variable reduces it.

::: proposition Toom–Cook with $N$ points {#prop-toom}
Let $N = \ell_1 + \ell_2 - 1$ and let $x_1, \ldots, x_N$ be distinct
points of $\Z \cup \{\infty\}$, where $a(\infty)$ denotes the leading
coefficient. Then $c = ab$ is determined by the $N$ products
$a(x_j) b(x_j)$, and it can be computed from them with additions,
multiplications and exact divisions by small integers.
:::

::: proof
$c$ has degree $N - 1$, so it is determined by its values at $N$
points (at $\infty$, by its leading coefficient), and $c(x_j) =
a(x_j) b(x_j)$. Solving the Vandermonde system over $\mathbb{Q}$ by
Newton's divided differences, or with a precomputed inverse, involves
divisions by differences of points; since the solution $c$ has integer
coefficients, the divisions can be arranged to be exact.
:::

With $N$ points the product needs $N$ multiplications instead of
$\ell_1 \ell_2$, which is optimal for bilinear algorithms. The
additions and the multiplications of the evaluation and interpolation
by powers of the points grow quadratically with $N$, so the method
pays off only when the coefficient products dominate. Toom used this
idea for integers [@Too1963], [@Coo1966]; Bodrato and Zanoni studied
optimal sequences of operations [@BZ2007]. `_fmpz_poly_mul_toom_scalar`
uses the points $0, \pm 1, \ldots, \pm p, \infty$ (and one more point
for odd $N$), up to $N = 20$ on 64-bit machines, with a table of
interpolation operations whose divisions are by single words.

Longer products split recursively: `_fmpz_poly_mul_toom_karatsuba`
applies Karatsuba's algorithm [@KO1962] to the polynomials, cutting the
longer operand into blocks of the length of the shorter for unbalanced
products, with Toom–Cook as base case for products of length up to 5,
11, 15 or 20 as the coefficients grow from below 2500 to over 6000
bits. Truncated products use Mulders' short product [@Mul2000] (high
products by reversal), and middle windows the Karatsuba middle product
of Hanrot, Quercia and Zimmermann [@HQZ2004], which computes the
middle product of a $2m \times m$ product in about the time of an $m
\times m$ product ([](algorithms/polynomials/middle-product)). The
Toom–Cook middle product counts multiplications in the same way: a
window of length $w$ of a product of operands of lengths $m, n$ needs
the sum of the two smallest of $m, n, w$, minus one.

## Kronecker substitution {#ks}

Kronecker substitution ([](algorithms/polynomials/kronecker)) packs
each operand into one integer by evaluating at $2^k$ and reads the
coefficients of $c$ from the integer product. Over $\Z$ the
coefficients have signs, which the packing represents by borrows.

::: lemma Signed digits {#lem-signed}
Let $k \ge 1$ and let $c = \sum_{j < N} c_j x^j \in \Z[x]$ with $-2^{k-1}
\le c_j < 2^{k-1}$ for all $j$. Then $c$ is determined by $C = c(2^k)$:
$c_0$ is the residue of $C$ modulo $2^k$ in $[-2^{k-1}, 2^{k-1})$, and
$(C - c_0)/2^k = (c - c_0)/x$ evaluated at $2^k$.
:::

::: proof
$C - c_0 = 2^k \sum_{j \ge 1} c_j 2^{k(j-1)}$, so $C \equiv c_0 \pmod
{2^k}$, and $c_0$ is the unique residue in the half-open interval of
length $2^k$. Induction on $N$ gives the other coefficients.
:::

`_fmpz_poly_mulmid_KS` takes $k = \beta + \lceil \log_2 (\ell_2 + 1)
\rceil + s$, where $s = 1$ if a coefficient is negative and 0
otherwise, so that $|c_j| < 2^{k - s}$ by the bound above. It packs
both operands with `_fmpz_poly_bit_pack`, negated if necessary so that
the packed integers are positive, multiplies with `flint_mpn_mul` (or
`flint_mpn_sqr`), and unpacks with `_fmpz_poly_bit_unpack`, which
propagates the borrows of [](#lem-signed); without negative
coefficients the fields are read as unsigned numbers. For a window it
computes only the needed limbs with `flint_mpn_mulmid`, with three
guard limbs below the window: a borrow that the truncated integer
product loses is absorbed by the slack in the field size, and the
unpacking corrects it. The cost is one integer product of about
$(\ell_1 + \ell_2) k$ bits, so KS is within a constant factor of the
lower bound whenever $\log \ell_2$ is small compared to $\beta$. The
integer product itself goes to `fft_small` above a few hundred limbs
([](algorithms/fft/choosing)).

## Multimodular FFT {#multimodular}

For long polynomials with small coefficients, KS wastes a large part of
each field on the $\log \ell_2$ bits of growth and on padding. The
small-prime FFT of `fft_small` instead reduces the coefficients modulo
2 to 8 primes of 50 bits, multiplies the images with number-theoretic
transforms, and reconstructs each coefficient by the Chinese remainder
theorem in the symmetric range, which is exact if the product $P$ of
the primes satisfies $P > 2 \ell_2 2^\beta$
([](algorithms/fft/multimodular)). With eight primes this covers $\beta
+ \log_2 \ell_2$ up to about 397 bits;
`_fmpz_poly_mul_mid_default_mpn_ctx` declines larger coefficients, and
the dispatcher then uses KS, whose integer product is computed by the
same transforms. See [](fft_small/poly_mul).

## Schönhage–Strassen convolutions {#ss}

`_fmpz_poly_mulmid_SS` maps each coefficient to an element of $R =
\Z/(2^N + 1)\Z$ with $N \ge \beta + \lceil \log_2 \ell_2 \rceil + 1$,
rounded up as the $\sqrt 2$ trick requires, and computes the product as
a convolution over $R$ of power-of-two length, in which every twiddle
factor is a power of 2 ([](algorithms/fft/schonhage-strassen)). Since
$|c_k| < 2^N/2$, the symmetric residues modulo $2^N + 1$ are the
coefficients. The pointwise products modulo $2^N + 1$ are themselves
computed by `fft_small` when a suitable ring exists
([](fft_small/mulmod_2expp1)). For windows that do not reach the top of
the product, a cyclic convolution of length $L \ge \max(n_{hi},
\ell_1 + \ell_2 - 1 - n_{lo})$ suffices, because the wrapped-around
terms fall below $n_{lo}$ ([](algorithms/fft/multimodular#thm-wrap));
`fft_convolution_choose_mid` uses it when it is shorter than the
truncated transform of the full product, which saves up to a third of
the work for the middle products of Newton iteration. The coefficient
ring grows with $\beta$ while the transform length depends only on the
lengths, which suits short polynomials with large coefficients, and the
convolutions parallelise well.

## Classical multiplication in transform space {#transform-space}

For short polynomials with very large coefficients, the integer
products dominate, and an $\ell_1 \times \ell_2$ classical product
multiplies each coefficient by $\ell_2$ others.
`_fmpz_poly_mulmid_classical_fft_small` transforms each coefficient
once with `fft_small` and computes each output coefficient as a sum of
pointwise products, transformed back once
([](fft_small/transform-space)). This needs one forward transform per
input coefficient and one inverse transform per output coefficient,
instead of three transforms per coefficient product, and it beats Karatsuba for truncated
products with coefficients of tens of thousands of bits up to a few
dozen terms.

## Selection {#selection}

The choice in `_fmpz_poly_mulmid` was tuned on x86-64 with `fft_small`
for random coefficients of uniform size. Roughly, with $n$ the length of
the shorter operand:

| | $\beta < 3000$ bits | $\beta \ge 3000$ bits |
|---|---|---|
| $n \le$ 2–15 | word accumulation or classical | classical (truncated products), Karatsuba–Toom |
| $n \lesssim$ 20–80 | packed limbs | Karatsuba–Toom up to 4–18 terms; transform-space classical for truncated products up to 21–87 terms when $\beta \ge 16\,000$ |
| larger | `fft_small` multimodular if $\beta \lesssim 400$, else KS | KS; SS for huge products or many threads |

Middle products with $n \ge 16$ and $\beta$ around 2000–4400 or above
120000 bits use SS when a short cyclic convolution exists, and operands
whose coefficient sizes differ by a factor of at least 8 use the
direct algorithms, which pay only for the actual sizes. The exact
conditions are listed in [](fmpz_poly/mul#selection). Builds without
`fft_small` use only KS and SS for large inputs, choosing by the
number of limbs per coefficient.

## Use in FLINT {#use}

- `_fmpz_poly_mulmid` dispatches; `_fmpz_poly_mul`, `_fmpz_poly_mullow`,
  `_fmpz_poly_mulhigh`, `_fmpz_poly_sqr`, `_fmpz_poly_sqrlow` and
  `fmpz_poly_mulhigh_n` are thin wrappers, as is the method
  `GR_METHOD_POLY_MULLOW` of the ring $\Z$ in `gr`.
- The algorithms are `_fmpz_poly_mulmid_classical`,
  `_fmpz_poly_mulmid_mpn`, `_fmpz_poly_mulmid_toom_scalar`,
  `_fmpz_poly_mulmid_toom_karatsuba`, `_fmpz_poly_mulmid_KS`,
  `_fmpz_poly_mul_mid_default_mpn_ctx` (in `fft_small`),
  `_fmpz_poly_mulmid_SS` and `_fmpz_poly_mulmid_classical_fft_small`.
- `fmpz_poly_mul_SS_precache_init` keeps the transform of a fixed
  operand for repeated products.
- Products in `fmpq_poly`, `fmpz_poly_q`, `fmpz_mod_poly`, `padic_poly`
  and `nf_elem` call `_fmpz_poly_mul` or `_fmpz_poly_mullow`, and
  `arb_poly` multiplies blocks of exact midpoints with them.

::: functions _fmpz_poly_mulmid _fmpz_poly_mulmid_classical _fmpz_poly_mulmid_mpn _fmpz_poly_mulmid_toom_karatsuba _fmpz_poly_mulmid_KS _fmpz_poly_mulmid_SS _fmpz_poly_mulmid_classical_fft_small :::
