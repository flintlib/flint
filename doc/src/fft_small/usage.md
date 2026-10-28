---
title: Where FLINT uses fft_small
summary: Which FLINT functions multiply through fft_small and from which sizes, how many primes a product needs, and what the default contexts cost in memory and threads.
order: 1
---

`fft_small` is FLINT's main engine for large products. Integer, polynomial
and matrix arithmetic reach it through their own dispatchers, which
compare the operand sizes with tuned thresholds; this guide collects those
thresholds and explains what a call costs. The module exists only when
`FLINT_HAVE_FFT_SMALL` is 1 in `flint-config.h`, which `configure` sets on
64-bit x86 with AVX2 and on 64-bit ARM with NEON. Elsewhere the same
functions multiply large integers with GMP's `mpn_mul` and, from
`FLINT_FFT_MUL_THRESHOLD` (32000) limbs, the Schönhage–Strassen code of
the [`fft`](fft) module.

## Callers and thresholds {#callers}

Thresholds depend on the target and are defined in `flint-mparam.h`. The
first value is that of the Intel targets from Skylake on, which an
installation for a recent x86-64 processor uses; the others are given in
parentheses.

| Function | Entry point | Used when |
|---|---|---|
| `flint_mpn_mul`, `flint_mpn_mul_n` | `mpn_mul_default_mpn_ctx` | $y_n \ge T$, or $x_n + y_n > 2T$ and $3 y_n > T$, for $x_n \ge y_n$; $T$ = `FLINT_FFT_SMALL_MUL_THRESHOLD` = 500 (240 on Zen 3 and 4, 810 on Apple M1, 1000 on other ARM, Broadwell and generic x86-64) |
| `flint_mpn_sqr` | `mpn_mul_default_mpn_ctx` | $n \ge$ `FLINT_FFT_SMALL_SQR_THRESHOLD` = 500 (400, 865, 1400) |
| `flint_mpn_mulhigh_n`, `flint_mpn_mullow_n` | `_mpn_ctx_mpn_mul_range` | $n >$ `FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF` = 500 (300 on Zen) |
| `flint_mpn_sqrhigh` | `_mpn_ctx_mpn_mul_range` | $n >$ `FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF` = 580 (470 on Zen) |
| `flint_mpn_mulmid` | `_mpn_ctx_mpn_mul_range` | $\min(a_n, b_n) \ge T$ (the squaring threshold if $a = b$), and windows of at least 48 limbs; a window shorter than the longer operand is computed as a transposed product when that has the shorter transform |
| `flint_mpn_mul_complex` and its variants | `_fft_small_mpn_mul_complex` and relatives (`_fft_small_mpn_mul_shared` for a real or imaginary operand times a complex one) | parts of at least 230 limbs (270 for squares) |
| `flint_mpn_mul_shared` | `_fft_small_mpn_mul_shared` | the longer product passes the test of `flint_mpn_mul` with $T$ = `FLINT_MPN_MUL_SHARED_THRESHOLD` and the shorter one with $T/2$; $T$ = `FLINT_FFT_SMALL_MUL_THRESHOLD` (240 on the `sapphirerapids` target) |
| `flint_mpn_gcd`, `flint_mpn_gcdext` (and through them `fmpz_gcd`, `fmpz_xgcd`, `fmpz_invmod`) | `fft_small_plan_init_mpn`, `fft_small_plan_init_mpn_cyclic` and the operand functions | the $2 \times 2$ matrix products of the half-GCD, from products in the range of `flint_mpn_mul` with $T$ replaced by $0.3\,T$ ([](algorithms/integers/gcd#use)) |
| `mp_real` binary splitting (hypergeometric series, `exp`, `sin`, `cos`) | `_fft_small_mpn_mul_shared`, `_fft_small_mpn_mul_precomp` | pairs of exact products with a common factor in the range of `flint_mpn_mul_shared`; powers of the argument used by at least two merges, for products in the FFT range |
| `nmod_poly_mul`, `_mullow`, `_mulmid`, `_sqr` | `_nmod_poly_mul_mid_default_mpn_ctx`; for middle products also `_nmod_poly_mul_tr_default_mpn_ctx` | lengths above a tuned table indexed by the bit size of $n$ (below); a middle product is computed as a transposed product when that has the shorter transform |
| `nmod_poly` evaluation, interpolation and extrapolation at geometric progressions | `_mul_precomp_init_tr` and `_nmod_poly_mul_tr_precomp` (evaluation, extrapolation), `_mul_precomp_init` and `_nmod_poly_mul_mid_precomp` (interpolation) | with stored transforms (bit 3 of `function`), where the corresponding product would use `fft_small` and the stored transform is cheaper |
| `fmpz_poly_mul`, `_mullow`, `_mulmid`, `_sqr` | `_fmpz_poly_mul_mid_default_mpn_ctx` | second length from 80 to 100 (100 to 130 for low products) depending on the coefficient sizes, unless they are very unequal; other methods take over if eight primes are too few |
| `fmpz_poly` and `fmpzi` Schönhage–Strassen | `sd_fft_mpn_mulmod_2expp1` | pointwise products in rings of at least 16384 bits |

The polynomial thresholds come from the table `fft_mul_tab` in
`nmod_poly/mullow_fft_small.c`: a full product of lengths $\ell_1 \ge
\ell_2$ uses `fft_small` when $\min(\ell_1, 2 \ell_2)$ reaches the entry
for the bit size of $n$, which is between 80 and 1326, for example 306 for
8 bits, 99 for 16 bits, 176 for 32 bits, 162 for 50 bits and 94 for 64
bits. The irregular entries reflect the Kronecker substitution variants
that compete with it. Squares and low products have their own tables;
moduli that are context primes use `fft_small` from $\min(\ell_1, 2
\ell_2) = 24$ (40 for low products), and moduli $n \le 23$ first pack two
coefficients into one (from total lengths of about 370). Middle products
use the threshold of the product they are the transpose of, and
`_nmod_poly_mulmid_fft_small` computes them as that
transposed product (`_nmod_poly_mul_tr_default_mpn_ctx`) when its
transform is shorter than that of the window of the full product
([](fft_small/poly_mul#transposed)).

Several other modules build their own drivers on the plans, operands and
the convolution engine of this module, or on the transforms directly:

- the transformed integer and polynomial rings of `gr`
  (`gr_ctx_init_transformed_mpn` and the transformed polynomial rings of
  `nmod` and `mpn_mod`), used by `fmpz_powm` and `fmpz_mod_pow` for
  moduli of at least 70000 bits, `fmpz_mat_mul_fft_small`,
  `nmod_poly_mat_mulmid_fft_small`, the matrix products of the
  polynomial half-gcd in `gr_poly`, and division with precomputed inverse
  over `nmod` from length 160;
- `mpn_mod` polynomial multiplication and the middle products of the
  decimal `radix` arithmetic, through `_fft_small_conv` (for `radix`,
  in its transposed mode for windows shorter than the longer operand,
  see `radix_mulmid_fft_small`);
- multipoint resultants in `nmod_mpoly` and `fmpz_mpoly_factor`, which
  transform with `sd_fft_trunc` modulo primes accepted by
  `fft_small_prime_supports_depth`, and `fmpz_poly_gcd_modular`, which
  tries the context primes first because products modulo them need a
  single transform.

## Which pipeline and how many primes {#primes}

For integers, operands of up to 40000 limbs use two primes and chunks of
40 to 45 bits; larger ones use four to eight primes and chunks of 84 to
192 bits ([](fft_small/mpn_mul)). For polynomials over $\Z/n\Z$ the
number of primes grows with the bit size of $n$ and, slowly, with the
length ([](fft_small/poly_mul#primes)): one prime for moduli of up to
about 16 bits, two up to about 44 bits, three above. Moduli that are
among the eight context primes, and primes $n$ of at least 20 bits with
a large power of two dividing $n - 1$, need a single transform without
Chinese remaindering. The cost of a product is roughly proportional to
the number of primes times the transform length.

The plan constructors make the same choices as the multiplication
drivers (for one thread) and can be used to inspect them:

```c example name=fft_small_choices requires=fft_small
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/fft_small.h>

int main(void)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    slong sizes[5] = { 1000, 40000, 100000, 1000000, 10000000 };
    int bits[4] = { 16, 32, 50, 64 };
    int i;

    /* integer products of two n-limb operands */
    for (i = 0; i < 5; i++)
    {
        fft_small_plan_t P;
        if (fft_small_plan_init_mpn(P, R, sizes[i], sizes[i], 1, 0))
            flint_printf("%8wd limbs: %wu primes, %3wu-bit chunks, length 2^%wu\n",
                         sizes[i], P->np, P->bits, P->depth);
        fft_small_plan_clear(P);
    }

    /* polynomials of length 10^5 modulo n = 2^k - 1 */
    for (i = 0; i < 4; i++)
    {
        fft_small_plan_t P;
        nmod_t mod;
        ulong len = 100000;

        nmod_init(&mod, (bits[i] == 64) ? UWORD_MAX : (UWORD(1) << bits[i]) - 1);
        if (fft_small_plan_init_nmod(P, R, 0, 2 * len - 1, 2 * len - 1,
                n_round_up(len, BLK_SZ), len, 2 * NMOD_BITS(mod), mod, 0))
            flint_printf("%2d-bit modulus: %wu prime%s, length 2^%wu\n",
                         bits[i], P->np, P->np == 1 ? "" : "s", P->depth);
        fft_small_plan_clear(P);
    }

    flint_cleanup_master();
    return 0;
}
```
```output
    1000 limbs: 2 primes,  44-bit chunks, length 2^12
   40000 limbs: 2 primes,  41-bit chunks, length 2^17
  100000 limbs: 5 primes, 116-bit chunks, length 2^17
 1000000 limbs: 6 primes, 136-bit chunks, length 2^20
10000000 limbs: 7 primes, 160-bit chunks, length 2^23
16-bit modulus: 1 prime, length 2^18
32-bit modulus: 2 primes, length 2^18
50-bit modulus: 3 primes, length 2^18
64-bit modulus: 3 primes, length 2^18
```

## Memory and threads {#resources}

Each thread that multiplies has its own default context, created on first
use: about 170 KB of tables, plus twiddle factors for transforms longer
than $2^{12}$, which take half a double per point and prime (4 MB per
prime at length $2^{20}$), plus a scratch buffer that grows to the
largest product computed so far and is kept for reuse. `flint_cleanup`
frees the default context of the calling thread, and
`flint_cleanup_master` also stops the thread pool; long-running programs
that once multiplied very large numbers can call them to return the
memory.

The multiplication functions use the threads of FLINT's pool, as set by
`flint_set_num_threads`. Integer products above 40000 limbs use up to
eight threads, polynomial products use one thread per prime from lengths
of about 1500 and up to eight threads for the reconstruction. The
two-prime integer pipeline, and therefore every product below 40000
limbs, is single-threaded. Threads are only borrowed for the duration of
a call, so the functions can be called from several threads at once,
each using its own default context.

## Calling the module directly {#direct}

Calling `fft_small` directly is useful to

- multiply without the size dispatch, for example to measure;
- keep the transform of an operand that multiplies many others
  ([](fft_small/precomp) for polynomials, [](fft_small/mpn_shared) for
  integers);
- compute a middle product as a transposed product
  ([](fft_small/poly_mul#transposed));
- add several products before transforming back, as in a dot product of
  polynomials or integers or a matrix product, which saves one inverse
  transform and one reconstruction per product
  ([](fft_small/transform-space));
- compute products modulo $2^N + 1$ ([](fft_small/mulmod_2expp1)) or with
  the transforms themselves ([](fft_small/sd_fft)).

Code that must also build without `fft_small` should test
`FLINT_HAVE_FFT_SMALL` and fall back to `flint_mpn_mul`,
`_nmod_poly_mul` or `_fmpz_poly_mul`.
