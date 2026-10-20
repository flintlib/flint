---
title: FFTs in FLINT
kind: topic
summary: FLINT has three FFT implementations — Schönhage–Strassen over $\Z/(2^N+1)\Z$, small-prime transforms in floating point, and generic DFTs over any ring. This page compares them, says which one each multiplication routine uses, and maps the techniques they share.
order: 0
---

A fast Fourier transform evaluates a polynomial of length $L$ at the $L$
powers of a root of unity of order $L$ in $O(L \log L)$ ring operations,
and the inverse transform interpolates. Multiplying two polynomials then
costs two transforms, $L$ pointwise products and one inverse transform,
and integers are multiplied the same way after cutting them into
coefficients. The ring in which the transform runs decides almost
everything else: how the roots of unity are represented, how large the
coefficients of the product may be, and how fast a butterfly is on
current hardware. FLINT has one implementation for each of three choices.

## The three implementations {#implementations}

| | [`fft`](fft) | [`fft_small`](fft_small) | [`gr_dft`](gr_dft) |
|---|---|---|---|
| Ring | $\Z/(2^N + 1)\Z$ | $\Z/p\Z$ for primes $p < 2^{50}$ with a large power of 2 dividing $p - 1$ | any `gr` ring with a suitable root of unity: $\Z/m\Z$, finite fields, complex balls, fixed point, exact numbers |
| Roots of unity | powers of $2$, and $\sqrt 2 = 2^{3N/4} - 2^{N/4}$, so every twiddle factor is a shift | precomputed tables of residues | precomputed tables in the ring |
| Arithmetic | multi-limb additions and shifts | integers stored in doubles, vector FMA ([](algorithms/simd/fp-modular-arithmetic)) | the ring's own operations |
| Lengths | powers of two, truncated | powers of two, truncated at block granularity | any $n \ge 1$ |
| Exact range | one ring element holds thousands of bits | 50 bits per prime; 2 to 8 primes with CRT | exact, or rigorous error bounds for complex balls |
| Main use | integer and `fmpz_poly` products in builds without `fft_small`; some large `fmpz_poly` convolutions | almost all large integer and polynomial products | `acb_dft`, Dirichlet character transforms, DFTs over general rings |
| Requires | nothing | AVX2 or NEON and 64-bit limbs | nothing |
| Algorithm pages | [](algorithms/fft/schonhage-strassen) | [](algorithms/fft/small-prime-fft), [](algorithms/fft/multimodular) | [](algorithms/fft/dft-algorithms), [](algorithms/fft/fixed-point-dft) |

The first two exist to multiply. The Schönhage–Strassen transform needs
no multiplications at all, only shifts and additions of multi-limb
numbers, and its ring grows with the operands, which gives the bound
$O(n \log n \log \log n)$ for $n$-bit products [@SS1971]. The small-prime
transform does the opposite: it keeps every number in the 53-bit
significand of a double, so that a butterfly is a handful of vector
instructions, and it recovers large coefficients by transforming modulo
several primes and combining the results with the Chinese remainder
theorem. On processors with wide FMA units the second approach wins at
every size FLINT is tuned for: the overview of the [`fft`](fft) module
measures `fft_small` at 1.5 to 2 times the speed of
`flint_mpn_mul_fft_main` where both apply. The third implementation is
not a multiplication engine but a general DFT library, for transforms of
any length over any ring and for numerical transforms with proven error
bounds.

## Which one a product uses {#dispatch}

Users rarely call a transform directly; the multiplication functions
choose by size.

- **Integers.** `flint_mpn_mul` uses FLINT's own small-size code and
  Karatsuba, then GMP's `mpn_mul`, and for large operands
  `mpn_mul_default_mpn_ctx` of `fft_small`, from a smaller operand of
  `FLINT_FFT_SMALL_MUL_THRESHOLD` limbs (500 on recent Intel processors;
  [](fft_small/usage#callers) lists the other targets). Builds without
  `fft_small` call `flint_mpn_mul_fft_main` of `fft` instead, from
  `FLINT_FFT_MUL_THRESHOLD` = 32000 limbs. Everything built on
  `flint_mpn_mul`, starting with `fmpz_mul`, inherits this choice.
- **Polynomials over $\Z/n\Z$.** `nmod_poly_mul` and its low, middle and
  high variants use classical multiplication, then Kronecker
  substitution, then `fft_small`, with cutoffs that depend on the bit
  size of $n$ ([](algorithms/polynomials/multiplication)). Kronecker
  substitution ends in an integer product, so its large cases also reach
  an FFT, through `flint_mpn_mul`.
- **Polynomials over $\Z$.** For most large inputs `fmpz_poly_mul`
  calls `_fmpz_poly_mul_mid_default_mpn_ctx`, a direct multimodular
  product with up to eight 50-bit primes of `fft_small`, when the
  coefficients of the product (plus the growth from the length) have up
  to about 397 bits, and Kronecker substitution, whose integer product
  also uses `fft_small`, beyond that. It uses Schönhage–Strassen
  convolutions (`fft_convolution`) for huge products, for some squarings
  and middle products, and with four or more threads
  ([](algorithms/polynomials/multiplication-integers)). The pointwise products of those convolutions are
  computed by `fft_small` again ([](fft_small/mulmod_2expp1)).
- **Matrices.** `fmpz_mat_mul_fft_small` and
  `nmod_poly_mat_mulmid_fft_small` transform each entry once and do the
  matrix products in transform space, so that $O(n^2)$ transforms serve
  $n^3$ pointwise products ([](algorithms/fft/multimodular#reuse)).
- **Complex numbers.** `acb_dft` and `acb_dirichlet_dft` are wrappers
  around `gr_dft`, which computes rigorous enclosures in fixed-point
  arithmetic when it can ([](gr_dft/acb)).

## Shared techniques {#techniques}

The implementations differ in their rings but share most of their
structure.

**Truncation.** A radix-2 transform for a product of length $t$ needs
length $L = 2^{\lceil \log_2 t \rceil}$, which wastes up to half of the
work just above a power of two. Van der Hoeven's truncated FFT computes
only the $t$ values needed and inverts from them in about
$(t/2) \log_2 L$ butterflies. Both multiplication engines use it.

::: summary algorithms/fft/truncated-fft :::

**Cache blocking.** A long transform is split as $L = L_1 L_2$ into
short transforms that fit in cache: Bailey's four-step or matrix Fourier
algorithm. `fft` uses it from a tuned size on, fused with the pointwise
products; `fft_small` applies the same decomposition recursively, in a
form without twiddle multiplications ([](algorithms/fft/small-prime-fft#blocking));
`gr_dft` offers it as one of its plan types.

::: summary algorithms/fft/matrix-fourier :::

**Several primes and the CRT.** When the coefficients of a product do
not fit modulo one prime, the product is computed modulo several and
reconstructed. Matrix multiplication over $\Z/n\Z$ uses the same idea
with floating-point matrix products
([](algorithms/linear-algebra/modular-matrix-kernels)).

::: summary algorithms/fft/multimodular :::

**Middle products and wrap-around.** A cyclic convolution of length $L$
gives the product modulo $x^L - 1$. When only a window of the product is
needed, the wrapped-around terms may be allowed to land outside the
window, which can save up to half of the transform length when a power
of two lies in the right place. The transposition principle gives the
saving for every length: `fft_small` has transposed versions of its
truncated transforms, with which a middle product costs a product of
the output length and the shorter operand. `nmod_poly_mulmid`,
`flint_mpn_mulmid` and `radix_mulmid` use them for windows shorter than
the longer operand when the transform is shorter
(`_nmod_poly_mul_tr_mpn_ctx`, `_mpn_ctx_mpn_mul_range`,
`radix_mulmid_fft_small`), and the evaluation and extrapolation at
geometric progressions keep the transposed transform of their fixed
polynomial (`_mul_precomp_init_tr`).

::: summary algorithms/polynomials/middle-product :::

**Kronecker substitution** turns polynomial products into integer
products, which is how the integer FFTs serve polynomial arithmetic.

::: summary algorithms/polynomials/kronecker :::

## Using the transforms directly {#direct}

Each module has a guide for those who need the transforms rather than a
finished product: [](fft/using) multiplies integers by hand with the
splitting, convolution and recombination functions of `fft`;
[](fft_small/transform-space) keeps operands transformed and accumulates
products in transform space with `fft_small`; and
[](gr_dft/getting-started) computes DFTs over $\Z/p\Z$ and the complex
numbers with `gr_dft`.
