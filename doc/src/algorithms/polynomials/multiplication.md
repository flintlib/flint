---
title: Multiplying polynomials over Z/nZ
kind: topic
summary: How `nmod_poly` multiplies: classical products built on modular dot products, Kronecker substitution (KS, KS2, KS4), and the small-prime FFT of `fft_small`, with the length and modulus-size cutoffs that choose between them for full, low, middle and high products.
order: 0
---

Almost every fast algorithm for polynomials over $\Z/n\Z$ (division,
GCD, evaluation, interpolation, composition, power series) spends its
time in polynomial multiplication. FLINT multiplies polynomials with
word-size modulus $n$ in one of three ways, chosen by the lengths of
the operands and the bit length $\beta$ of $n$ (`NMOD_BITS`):

| Method | Entry points | Cost | Best for |
|---|---|---|---|
| Classical | `_nmod_poly_mul_classical` | $\ell_1 \ell_2$ word products, $\ell_1 + \ell_2 - 1$ reductions | short operands |
| [Kronecker substitution](algorithms/polynomials/kronecker) | `_nmod_poly_mul_KS`, `_nmod_poly_mul_KS2`, `_nmod_poly_mul_KS4` | one, two or four integer products, $(\ell_1 + \ell_2)(2\beta + \lg \ell_2)$ bits in total | medium lengths |
| Small-prime FFT | `_nmod_poly_mullow_fft_small`, `_nmod_poly_mulmid_fft_small` | $O(N \log N)$ for $N = \ell_1 + \ell_2$, with one to four primes | long operands |

Here $\ell_1 \ge \ell_2$ are the operand lengths. The algorithms are
standard; see [@GG2003, Ch. 8] and [@Ber2008] for surveys.

Products come in four shapes, each with its own selection among these
methods:

- the full product $fg$ (`_nmod_poly_mul`);
- the low product $fg \bmod x^m$ (`_nmod_poly_mullow`), used by Newton
  iteration and power series arithmetic;
- the middle product, coefficients $m_1 \le k < m_2$ of $fg$
  (`_nmod_poly_mulmid`), the transpose of multiplication, see
  [](algorithms/polynomials/middle-product);
- the high product, coefficients $k \ge m$ (`_nmod_poly_mulhigh`).

## Generic version {#generic}

`gr_poly_mul`, `gr_poly_mullow` and `gr_poly_mulmid` work over any ring
by calling the ring methods `GR_METHOD_POLY_MULLOW` and
`GR_METHOD_POLY_MULMID`. The `gr` context for $\Z/n\Z$
(`gr_ctx_init_nmod`) installs `_nmod_poly_mullow` and
`_nmod_poly_mulmid` as these methods, so generic code over this ring
(Newton iteration, subproduct trees, half-gcd, composition) makes
exactly the choices described here. Other rings have their own
methods: `fmpz_poly` products for $\Z$ and $\mathbb{Q}$, Kronecker substitution
for finite fields $\mathbb{F}_q$, block multiplication for real and
complex balls, multi-word kernels and the small-prime FFT for
`mpn_mod` ([](mpn_mod/poly-mul)); the table in
[](algorithms/polynomials/overview#multiplication) lists them. A ring
without a method of its own gets the generic default, the classical
algorithm with one dot product (`_gr_vec_dot_rev`) per output
coefficient, which is the analogue of the classical method below with
the ring's dot product in place of the modular one. `gr_poly` also
offers one level of Karatsuba, Toom-3 and Toom-Cook multiplication and
bivariate Kronecker substitution ([](gr_poly#mul-algorithms)) for rings
that want to install them.

## Classical multiplication {#classical}

The schoolbook method computes each output coefficient
$h_k = \sum_{i+j=k} f_i g_j$ as a dot product of a slice of $f$ with a
reversed slice of $g$, using `_nmod_vec_dot_rev`. As explained in
[](algorithms/linear-algebra/dot-products), the sum is accumulated
exactly in one, two or three words and reduced once, so a product of
lengths $\ell_1, \ell_2$ costs $\ell_1 \ell_2$ word multiplications but
only $\ell_1 + \ell_2 - 1$ modular reductions. When a whole product
coefficient fits in one word ($\ell_2 (n-1)^2 < 2^{\bits}$),
`_nmod_poly_mul_classical` instead accumulates all coefficients with
plain word arithmetic, row by row, and reduces at the end; squares are
computed with about half the multiplications. Classical low, middle and
high products (`_nmod_poly_mullow_classical`,
`_nmod_poly_mulmid_classical`, `_nmod_poly_mulhigh_classical`) compute
only the requested dot products, so they cost proportionally less.

FLINT does not use Karatsuba's method [@KO1962] for `nmod_poly`: Kronecker
substitution reaches the Karatsuba and Toom–Cook range of integer
multiplication with less overhead.

## Kronecker substitution {#ks}

Packing the coefficients into integers with fields of
$2\beta + \lg \ell_2$ bits turns the product into one integer
multiplication, which `flint_mpn_mul` performs with basecase, Toom–Cook
or FFT code as appropriate. Harvey's variants KS2 and KS4 evaluate at
$\pm 2^b$ (and at $\pm 2^{-b}$) and replace this product by two or four
products of a half or a quarter of the size. Details and proofs are in
[](algorithms/polynomials/kronecker).

## The small-prime FFT {#fft-small}

For long operands, `nmod_poly` uses the `fft_small` module, which
computes convolutions modulo primes $p < 2^{50}$ with truncated FFTs in
double-precision floating point, vectorised with AVX2, AVX-512 or NEON
([](algorithms/fft/small-prime-fft), [](algorithms/fft/truncated-fft)).
The coefficients $h_k \le \ell_2 (n - 1)^2$ of the integer product are
reconstructed by Chinese remaindering from their residues modulo
$s$ primes ([](algorithms/fft/multimodular)) and then reduced modulo $n$.
The default primes $p_1, \ldots, p_8$ have about 50 bits each; their
products have $49.98$, $99.21$, $149.15$ and $198.85$ bits for
$s = 1, 2, 3, 4$. The number of primes is the least $s$ with
$p_1 \cdots p_s \ge \ell_2 \, 2^{2\beta}$, so

| Modulus | Primes for typical lengths |
|---|---|
| $2\beta + \log_2 \ell_2 \le 49.98$ (e.g. $\beta \le 20$ and $\ell_2 \le 512$) | 1 |
| $2\beta + \log_2 \ell_2 \le 99.21$ (e.g. $\beta \le 40$ and $\ell_2 \le 2^{19}$) | 2 |
| $2\beta + \log_2 \ell_2 \le 149.15$ (e.g. $\beta = 64$, $\ell_2 \le 2^{21}$) | 3 |
| larger | 4 |

Three cases avoid Chinese remaindering altogether:

- If $n$ is itself one of the FFT primes of the default context, the
  transform is done modulo $n$ directly. For such $n$ (50 bits),
  `_nmod_poly_mullow_want_fft_small` selects the FFT from length 24
  instead of using the table below.
- If $n$ is any prime with $2^{19} \le n \le 2^{50}$, $2^d \mid n - 1$
  for the transform length $2^d$, $\ell_2 \ge 100$, and $n$ passes the
  error bound check of the floating-point arithmetic, the transform is
  done modulo $n$ with a temporary context (the primality test is cached
  per thread).
- If $n \le 23$, two coefficients are packed into each FFT slot, which
  halves the transform length (see [](#repack)).

The products are computed by `_nmod_poly_mul_mid_default_mpn_ctx` in
`fft_small` ([](fft_small/poly_mul)), which handles full, low and
middle products: for a middle product it chooses a cyclic convolution
length so that the wrap-around only pollutes coefficients outside the
requested range when a power of two allows it
([](algorithms/polynomials/middle-product#fft)). A middle product can
also be computed as a transposed product, at the cost of a product of
the shorter operand and the output length, with the transposed
transforms of `fft_small` (`_nmod_poly_mul_tr_default_mpn_ctx`);
`_nmod_poly_mulmid_fft_small` takes whichever of the two has the shorter
transform ([](algorithms/polynomials/middle-product#transposed-fft)).
If FLINT's thread pool
has workers and the shorter operand has length at least 1500, the
transforms for the different primes run in parallel threads (and from
length 5000 the two operands are transformed in parallel).
[](fft_small/usage) lists the other callers of `fft_small` and their
thresholds, and how to call the module directly.

### Repacking for tiny moduli {#repack}

For $n \le 23$ a product coefficient has few bits, and a 50-bit prime
has room for three base-$M$ digits with $M = 2^{16}$ or $2^{17}$. Write
$A(x) = A_e(x^2) + x A_o(x^2)$ and $B(x) = B_e(x^2) + x B_o(x^2)$ and
set $A_2 = A_e + M A_o$, $B_2 = B_e + M B_o$, polynomials of half the
length with coefficients below $M^2$.

::: lemma {#lem-repack}
$A_2 B_2 = A_e B_e + M (A_e B_o + A_o B_e) + M^2 A_o B_o$, and
$AB = (A_e B_e)(x^2) + x^2 (A_o B_o)(x^2) + x (A_e B_o + A_o B_e)(x^2)$.
If all coefficients of $A_e B_e$ and of $A_e B_o + A_o B_e$ are less than
$M$ and all coefficients of $A_2 B_2$ are less than the prime $p$, then
coefficient $k$ of $A_2 B_2$, computed modulo $p$, has base-$M$ digits
$u_k$, $v_k$, $w_k$ with $(AB)_{2k} = u_k + w_{k-1}$ and
$(AB)_{2k+1} = v_k$.
:::

The proof is by expanding the products. `_nmod_poly_mulmid_fft_small_repack`
checks the hypotheses with the bound $\min(\ell_1, \ell_2) (n-1)^2 < 2^{16}$
or, failing that, with Cauchy–Schwarz bounds
$\|A_e\|_2 \|B_e\|_2 < M$ and so on, which are much smaller for random
inputs. If neither works it returns 0 and the ordinary FFT is used.

::: note
In FLINT 3.6 and earlier,
the Cauchy–Schwarz test added the last coefficient of an operand of odd
length to $\|A_e\|_2^2$ (respectively $\|B_e\|_2^2$) without squaring
it, so the norm could be underestimated and $M = 2^{17}$ accepted
although a coefficient of $A_e B_e$ reaches $2^{17}$; with $n = 23$ and
$a = b = 22 \sum_{i=0}^{270} x^{2i}$, `nmod_poly_mul` was wrong. A check
is given in [](nmod_poly/mul#small-moduli).
:::

## Cutoffs {#cutoffs}

The following describes x86-64 with `fft_small` (FLINT built with AVX2
or AVX-512); without `fft_small` the FFT branch is absent, see below.
"Length" means $\min(\ell_1, 2\ell_2)$ for full products, so that very
unbalanced products are treated like balanced ones of length $2\ell_2$.

**Full products.** `_nmod_poly_mul`, with $\ell_1 \ge \ell_2$, tries
in order:

1. classical if $\ell_2 \le 5$;
2. the FFT if `_nmod_poly_mullow_want_fft_small` says so: for $n \le 23$
   when $\ell_2 \ge 64$ and $370 < \ell_1 + \ell_2 - 1$, except for
   lengths $513$ to $659$ (where the transform length doubles) and, for
   $n \le 3$, $1025$ to $1099$; for other $n$ when the length is at least
   a tabulated value for $\beta$ (separate tables for products and
   squares), or at least 24 for the 50-bit FFT primes;
3. classical if $3 \cdot \text{length} < 2 \max(\beta, 10)$;
4. KS if $\text{length} \cdot \beta < 800$;
5. KS2 if $\text{length} \cdot (\beta + 1)^2 < 100000$;
6. KS4 otherwise.

For balanced operands of length $\ell$ this gives the following ranges:

| $\beta$ | classical | KS | KS2 | KS4 | FFT (product) | FFT (square) | FFT primes at the cutoff |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 ($n \ge 24$) | $\le 6$ | 7–159 | 160–673 | – | $\ge 674$ | $\ge 689$ | 1 |
| 8 | $\le 6$ | 7–99 | 100–305 | – | $\ge 306$ | $\ge 353$ | 1 |
| 12 | $\le 7$ | 8–66 | 67–191 | – | $\ge 192$ | $\ge 279$ | 1 |
| 16 | $\le 10$ | 11–49 | 50–98 | – | $\ge 99$ | $\ge 159$ | 1 |
| 20 | $\le 13$ | 14–39 | 40–81 | – | $\ge 82$ | $\ge 139$ | 1 |
| 22 | $\le 14$ | 15–36 | 37–189 | 190–437 | $\ge 438$ | $\ge 626$ | 2 |
| 24 | $\le 15$ | 16–33 | 34–159 | 160–323 | $\ge 324$ | $\ge 569$ | 2 |
| 28 | $\le 18$ | 19–28 | 29–118 | 119–267 | $\ge 268$ | $\ge 292$ | 2 |
| 32 | $\le 21$ | 22–24 | 25–91 | 92–175 | $\ge 176$ | $\ge 191$ | 2 |
| 40 | $\le 26$ | – | 27–59 | 60–91 | $\ge 92$ | $\ge 159$ | 2 |
| 46 | $\le 30$ | – | 31–45 | 46–80 | $\ge 81$ | $\ge 200$ | 2 |
| 47 | $\le 31$ | – | 32–43 | 44–176 | $\ge 177$ | $\ge 191$ | 3 |
| 50 | $\le 33$ | – | 34–38 | 39–161 | $\ge 162$ | $\ge 182$ | 3 |
| 56 | $\le 37$ | – | – | 38–163 | $\ge 164$ | $\ge 182$ | 3 |
| 60 | $\le 39$ | – | – | 40–164 | $\ge 165$ | $\ge 182$ | 3 |
| 61 | $\le 40$ | – | – | 41–94 | $\ge 95$ | $\ge 152$ | 3 |
| 64 | $\le 42$ | – | – | 43–93 | $\ge 94$ | $\ge 145$ | 3 |

For squares the KS ranges extend up to the square cutoff. The jumps in
the FFT cutoff follow the number of primes (one prime up to
$\beta = 21$, two up to $46$, three above) and the word sizes of the KS
fields: for example, the drop at $\beta = 61$ coincides with the point
where the $2b$-bit digits of KS4 no longer fit in one word at these
lengths. The tables are compiled into the library; unlike some other
cutoffs in FLINT they do not depend on the processor.

**Low products.** `_nmod_poly_mullow` truncates the operands to the
output length $m$, calls `_nmod_poly_mul` if $m = \ell_1 + \ell_2 - 1$,
and otherwise uses classical for $\ell_2 \le 5$, the FFT if
$\min(\ell_1, \ell_2)$ reaches a third table (for $n \le 23$: $\ell_2 \ge 64$
and $\ell_1 + \ell_2 - 1 > 450$ except $513$ to $699$; for the 50-bit FFT
primes: $\min(\ell_1, 2\ell_2) \ge 40$), classical if
$m < 10 + \beta^2/10$, and `_nmod_poly_mullow_KS` otherwise. The tabulated
FFT cutoffs for low products decrease from 1115 for $\beta = 1$ to 69 for
$\beta = 21$, and lie between 174 and 220 for $\beta \ge 46$. Since KS
computes the whole product and the FFT needs a transform as long as for
the full product, a truncated product is substantially cheaper than a
full one only with the classical method. Low products use only plain
KS, not KS2 or KS4, which is one reason why the FFT cutoffs
differ from those for full products.

**Middle products.** `_nmod_poly_mulmid` first drops the operand
coefficients that do not contribute, puts the longer operand first, and
then uses classical if either $\ell_2$ or the output length $m$ is at
most 5, and `_nmod_poly_mul` if the full product is requested. Otherwise
it uses the FFT if `_nmod_poly_mulmid_want_fft_small` holds: if the FFT
would be used for a full product of lengths $\ell_2$ and $m$ (a middle
product is the transpose of such a product, and is computed as one when
that is shorter than the window of the full product), and also if
$\min(\max(\ell_2, m), 2\min(\ell_2, m)) \ge 100$ for
$9 \le \beta \le 15$ and $\beta \ge 32$. Below that it uses `_nmod_poly_mulmid_KS`
for small moduli and long outputs ($\ell_2 \ge 8$, and $m \ge \ell_2$ with
$\beta \le 8$ or $m \ge 8 \ell_2$ with $\beta \le 16$) and the classical
method otherwise.

**High products.** `_nmod_poly_mulhigh` uses the classical method for
tiny inputs ($\ell_1 + \ell_2 \le 6$, or $\ell_1 + \ell_2 < 16$ when a
coefficient fits in a word) and otherwise `_nmod_poly_mulmid` with the
window $[m, \ell_1 + \ell_2 - 1)$, which chooses its method as above, and
sets the coefficients below $m$ to zero.

**Without fft_small.** On machines without AVX2 or NEON support,
`_nmod_poly_mul` uses classical, KS, KS2 and KS4 by steps 1 and 3 to 6
above, so KS4 is used for all long products. Its integer products use
`flint_mpn_mul`, that is, GMP's `mpn_mul` (with its Toom–Cook code)
above the basecase and Karatsuba sizes and, above
`FLINT_FFT_MUL_THRESHOLD` limbs, FLINT's
Schönhage–Strassen FFT [@SS1971] ([](algorithms/fft/schonhage-strassen)).
Low and middle products use classical below $10 + \beta^2/10$ output
coefficients and KS above.

## Which function to call {#advice}

- Call `nmod_poly_mul`, `nmod_poly_mullow` or `nmod_poly_mulmid` (or
  their underscore forms); they make the choices above. The forced
  variants (`nmod_poly_mul_classical`, `nmod_poly_mul_KS`,
  `nmod_poly_mul_KS2`, `nmod_poly_mul_KS4`) are meant for testing and
  tuning.
- Ask only for the coefficients you need. A low product saves about half
  of the work with the classical method (with the FFT and KS it costs
  about as much as the full product). A middle product of lengths
  $m + n - 1$ and $n$ costs about as much as a product of lengths $m$ and
  $n$, with every method except KS; Newton iteration in FLINT relies on
  this ([](algorithms/polynomials/newton)).
- To multiply many polynomials by the same polynomial modulo the same
  modulus, precompute: `nmod_poly_mulmod_precond` stores data for a fixed
  factor, and the geometric progression structures of
  [](algorithms/polynomials/subproduct-tree#geometric) store FFTs of the
  fixed factors.
- Squares are detected when both operands are the same pointer with the
  same length; they use separate cutoffs and are cheaper by a factor of
  up to two.

## Other modules {#other}

`fmpz_poly` multiplies integer polynomials with a much larger set of
methods (classical, Karatsuba, Toom, KS, Schönhage–Strassen and the
multimodular `fft_small` product), all selected in `_fmpz_poly_mulmid`;
see [](algorithms/polynomials/multiplication-integers).
Polynomials over finite fields $\mathbb{F}_q$ are mostly multiplied by
Kronecker substitution into $\Z[x]$ or $\mathbb{F}_p[x]$
([](algorithms/polynomials/kronecker#use)), so they too end up in the
code described here.

::: functions _nmod_poly_mul _nmod_poly_mullow _nmod_poly_mulmid _nmod_poly_mulhigh _nmod_poly_mullow_want_fft_small _nmod_poly_mulmid_want_fft_small :::
