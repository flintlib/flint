---
title: Multiplying integers and polynomials
summary: When to call the module at all, how to multiply integers with it, and how to set up the coefficient arrays, ring and scratch space for convolutions of integer polynomials.
order: 1
---

The module offers three levels of functionality: integer multiplication
with automatically chosen parameters (`flint_mpn_mul_fft_main`),
convolutions of vectors of large integers with parameters chosen by the
caller (`fft_convolution` and its variants), and the transforms
themselves.

## Choosing a function {#choosing}

- To multiply integers, call `flint_mpn_mul` or `fmpz_mul`. For large
  operands they use `mpn_mul_default_mpn_ctx` in builds with `fft_small`
  (on x86-64 with AVX2 and on ARM64), and otherwise
  `flint_mpn_mul_fft_main`. Call `flint_mpn_mul_fft_main` directly only
  to get this algorithm regardless of the build, for example for
  comparisons; see [](fft/mul).
- To multiply polynomials with integer coefficients, call `fmpz_poly_mul`
  and its relatives, which choose between Schönhage–Strassen convolutions
  from this module, Kronecker substitution and classical methods. For many
  products with the same operand, `fmpz_poly_mul_SS_precache_init` and
  `fmpz_poly_mullow_SS_precache` precompute its transform.
- Call `fft_convolution` directly to multiply vectors of large integers
  that are not `fmpz_poly` objects, to control the transform parameters,
  or to multiply modulo $2^N + 1$ coefficientwise ([](fft/convolution),
  [](fft/precache)).
- For a single product modulo $2^N + 1$, call `fft_mulmod_2expp1` with a
  limb count from `fft_adjust_limbs` ([](fft/mulmod)).

## Memory layout {#layout}

A convolution of transform length $4n$, $n = 2^{\mathit{depth}}$, with
coefficients modulo $p = 2^N + 1$, $N = \mathit{limbs} \cdot \bits$,
needs:

1. For each operand an array of $4n$ pointers, each to its own block of
   $\mathit{limbs} + 1$ limbs. One coefficient is a signed two's
   complement integer in this block.
2. For each thread (`flint_get_num_threads()` of them), three blocks of
   $\mathit{limbs} + 1$ limbs, collected in the pointer arrays `t1`, `t2`
   and `s1`, and one block of $2(\mathit{limbs} + 1)$ limbs in the array
   `tt`.

The transforms swap pointers between the operand arrays and the scratch
arrays, so all blocks should come from a single allocation, freed as a
whole at the end, and the coefficients must always be accessed through
the pointer array.

## Choosing the parameters {#parameters}

For a product with at most $\ell$ coefficients whose absolute values are
below $2^{N - 1}$ (or nonnegative and below $2^N$):

1. Take $\mathit{depth} = \lceil \log_2 \ell \rceil - 2$ so that
   $2n < \ell \le 4n$, and pass $\mathit{trunc} = \ell$. Products of at
   most 2 coefficients are not supported.
2. Round $N$ up to a multiple of $n$ and of $\bits$, so that every point
   of the transform carries $w = N/n$ bits.
3. If $\mathit{limbs} = N/\bits$ exceeds 128, replace it by
   `fft_adjust_limbs(limbs)`, or by `fft_adjust_limbs_sd_fft(limbs, n, &m)`
   to allow the faster pointwise engine of `fft_small`.

The coefficient bound comes from the operands: with input coefficients
below $2^{b_1}$ and $2^{b_2}$ in absolute value and the shorter operand
of length $\ell_2$, the product coefficients are below
$\ell_2 \, 2^{b_1 + b_2}$. When the inputs are chunks of integers, as
below, the chunk size $b$ is chosen instead so that
$2b + \mathit{depth} + 1 \le N$; this is the rule of
`mul_truncate_sqrt2`.

## A complete integer product {#example}

The following program multiplies two integers of about 6300 bits with a
convolution of length 256 modulo $2^{128} + 1$, which is what
`mul_truncate_sqrt2` does with $\mathit{depth} = 6$ and $w = 2$. The
chunks have $b = \lfloor (128 - 7)/2 \rfloor = 60$ bits, so the product
coefficients are below $105 \cdot 2^{120} < 2^{128}$.

```c example name=fft_pipeline
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/mpn_extras.h>

int main(void)
{
    /* transform length 4n = 256, coefficients modulo 2^128 + 1 */
    slong depth = 6, n = WORD(1) << depth, limbs = 2, size = limbs + 1;
    flint_bitcnt_t bits = (limbs * FLINT_BITS - depth - 1) / 2;   /* 60 */
    slong nt = flint_get_num_threads(), i, j1, j2, trunc;
    mp_limb_t ** ii, ** jj, ** t1, ** t2, ** s1, ** tt, * mem, * r;
    mp_size_t n1, n2;
    mpz_t a, b, c;

    mpz_init(a); mpz_init(b); mpz_init(c);
    mpz_ui_pow_ui(a, 3, 4000);
    mpz_ui_pow_ui(b, 5, 2700);
    n1 = mpz_size(a);
    n2 = mpz_size(b);

    /* 2 x 4n coefficients, plus per-thread scratch: t1, t2, s1 of
       size limbs + 1 and tt of size 2 (limbs + 1) */
    ii = flint_malloc(8 * n * sizeof(mp_limb_t *));
    jj = ii + 4 * n;
    t1 = flint_malloc(4 * nt * sizeof(mp_limb_t *));
    t2 = t1 + nt; s1 = t2 + nt; tt = s1 + nt;
    mem = flint_malloc((8 * n + 5 * nt) * size * sizeof(mp_limb_t));
    for (i = 0; i < 8 * n; i++)
        ii[i] = mem + i * size;
    for (i = 0; i < nt; i++)
    {
        t1[i] = mem + (8 * n + i) * size;
        t2[i] = t1[i] + nt * size;
        s1[i] = t2[i] + nt * size;
        tt[i] = mem + (8 * n + 3 * nt + 2 * i) * size;
    }

    /* split into 60-bit coefficients and pad with zeros */
    j1 = fft_split_bits(ii, mpz_limbs_read(a), n1, bits, limbs);
    j2 = fft_split_bits(jj, mpz_limbs_read(b), n2, bits, limbs);
    for (i = j1; i < 4 * n; i++) flint_mpn_zero(ii[i], size);
    for (i = j2; i < 4 * n; i++) flint_mpn_zero(jj[i], size);

    trunc = j1 + j2 - 1;            /* must satisfy 2n < trunc <= 4n */
    flint_printf("j1 = %wd, j2 = %wd, trunc = %wd, 4n = %wd\n", j1, j2, trunc, 4 * n);

    fft_convolution(ii, jj, depth, limbs, trunc, t1, t2, s1, tt);

    /* evaluate the product polynomial at 2^bits */
    r = flint_calloc(n1 + n2, sizeof(mp_limb_t));
    fft_combine_bits(r, ii, trunc, bits, limbs, n1 + n2);

    mpz_mul(c, a, b);
    flint_printf("equal to mpz_mul: %d\n",
        mpn_cmp(r, mpz_limbs_read(c), mpz_size(c)) == 0 &&
        (mpz_size(c) == (size_t) (n1 + n2) || r[n1 + n2 - 1] == 0));

    flint_free(r); flint_free(mem); flint_free(t1); flint_free(ii);
    mpz_clear(a); mpz_clear(b); mpz_clear(c);
    return 0;
}
```
```output
j1 = 107, j2 = 105, trunc = 211, 4n = 256
equal to mpz_mul: 1
```

The same program works for larger transforms, where `fft_convolution`
switches to the threaded matrix Fourier algorithm; with
$\mathit{depth} = 9$ and eight limbs it multiplies integers of half a
million bits. The scratch arrays must then really have one entry per
thread: the number of threads must not be increased between allocating
them and the call.

## Polynomials with integer coefficients {#polynomials}

For polynomials, each coefficient goes into its own block instead of
being cut from an integer. Negative coefficients are stored in two's
complement over all $\mathit{limbs} + 1$ limbs, and a coefficient of the
result is negative if its residue exceeds $p/2$; `_fmpz_vec_get_fft` and
`_fmpz_vec_set_fft` convert between `fmpz` vectors and this format. The
example in [](fft/convolution#examples) multiplies two polynomials with
signed coefficients this way and compares the result with
`fmpz_poly_mul`; the second example there uses the pointwise engine of
`fft_small`. When only some coefficients of the product are wanted,
`fft_convolution_choose_mid` can select a shorter cyclic transform.

## Threads {#threads}

`fft_convolution` with $\mathit{depth} > 6$, `fft_convolution_sd_fft`,
`mul_mfa_truncate_sqrt2` (and so `flint_mpn_mul_fft_main` for large
operands) and the splitting functions run in parallel on the threads of
FLINT's global thread pool, as set by `flint_set_num_threads`. The
precached convolutions without the `fft_small` engine and, apart from
the splitting, `mul_truncate_sqrt2` run in the calling thread.
