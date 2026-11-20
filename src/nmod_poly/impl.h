/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_IMPL_H
#define NMOD_POLY_IMPL_H

#include "nmod_types.h"
#include "ulong_extras.h"  /* n_mulhi */
#include "nmod_vec.h"
#include "flint-mparam.h"  /* FLINT_NMOD_POLY_IFMA_* */

/*
    Products with a short operand by rows (_nmod_poly_mulmid_rows): one
    _nmod_vec_scalar_addmul_nmod of the longer operand per coefficient of
    the shorter one, of length len2 from 3 to _NMOD_POLY_ROWS_MAX_LEN2,
    rather than one dot product of at most len2 terms per output coefficient
    as in the classical code. _nmod_poly_mul, _nmod_poly_mullow and
    _nmod_poly_mulmid use them in place of the classical code and of KS when
    these scalar products have a vector kernel for the modulus and every row
    is long enough: the shortest row, of length row_len, reaches the length
    of the kernels and ratio * len2, since a row costs a few times more than
    a dot product on top of their common arithmetic. For powers of 2, whose
    dot products need no reduction, only up to 2^31 (2^48 with the 512-bit
    kernels) and with rows twice as long, 2 ratio len2: there the rows are
    1.1-4x faster than the classical code, which wins on shorter rows and
    from 2^32 with AVX2 (the fp50 kernel) or 2^50 with AVX-512 (measured Nov
    2026).

    _NMOD_POLY_ROWS_RATIO: the ratio for full and low products, 4 with the
    AVX2 kernels above 32 bits (doubles, peeling to an aligned address,
    scalar tails), 2 with the 512-bit kernels (the only ones above 50 bits,
    masked tails) and up to 32 bits.

    _NMOD_POLY_ROWS_MAX_LEN2: the largest len2, 12 with the NEON kernels
    (on the M4, at 31-32 bits, where the classical code does not
    accumulate in one word, the rows are 1.2-1.6x faster up to len2 =
    12-13, measured Nov 2026).

    _NMOD_POLY_NEON: 1 with the NEON kernels (aarch64), where KS compares
    better with the classical code in one word (on the M4, see mul.c,
    mullow.c and _nmod_poly_mulmid_want_one_word).
*/
#define _NMOD_POLY_ROWS_RATIO(mod) \
    ((NMOD_VEC_SIMD_SCALAR_MAX_BITS(1) > 50 || NMOD_BITS(mod) <= 32) ? 2 : 4)

#if defined(__ARM_NEON) && defined(__aarch64__)
# define _NMOD_POLY_NEON 1
# define _NMOD_POLY_ROWS_MAX_LEN2 12
#else
# define _NMOD_POLY_NEON 0
# define _NMOD_POLY_ROWS_MAX_LEN2 16
#endif

static inline int
_nmod_poly_want_rows(slong row_len, slong len2, slong ratio, nmod_t mod)
{
    return len2 >= 3 && len2 <= _NMOD_POLY_ROWS_MAX_LEN2
        && NMOD_BITS(mod) <= NMOD_VEC_SIMD_SCALAR_MAX_BITS(1)
        && row_len >= FLINT_MAX(NMOD_VEC_SIMD_MIN_LEN, ratio * len2)
        && ((mod.n & (mod.n - 1)) != 0
            || (NMOD_BITS(mod) <= (NMOD_VEC_SIMD_SCALAR_MAX_BITS(1) > 50 ? 49 : 32)
                && row_len >= 2 * ratio * len2));
}

void _nmod_poly_mulmid_rows(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

/* whether sums of len2 products of coefficients fit in one word: (n - 1)^2
   len2 < 2^64, which implies n <= 2^32 for len2 >= 2 */
static inline int
_nmod_poly_one_word(slong len2, nmod_t mod)
{
    return NMOD_BITS(mod) <= FLINT_BITS / 2
        && n_mulhi((mod.n - 1) * (mod.n - 1), len2) == 0;
}

/*
    Without NEON, the classical code accumulating in one word beats KS (or
    KS2, KS4) for a shorter operand of len2 < c(b) coefficients, b the
    number of bits of n, when the longer one has at least 2 len2 (full
    products, _nmod_poly_mul) or is cut to n >= 2 len2 (low products): c(b)
    = max(12, 2b) below 22 bits, 4b - 24 from 22 bits (fft_small takes
    over at cutoff_len 80-100 up to 21 bits, about 300 above). Measured
    crossovers on a Cascade Lake, gcc / clang, len1 >= 4 len2 (Nov 2026):
    12 / 9 at 4 bits, 16 / 12 at 8 bits, 18-32 / 18-20 at 12 bits, up to
    fft_small at 18-20 bits with gcc, 72-112 / 36-96 at 22-28 bits; Ice
    Lake and Zen 4 (without their IFMA kernel) and low products n x n/4 at
    least as far.
*/
static inline slong
_nmod_poly_one_word_cutoff(slong bits)
{
    return (bits < 22) ? FLINT_MAX(12, 2 * bits) : 4 * bits - 24;
}

/*
    Middle products by rows accumulating in one word
    (_nmod_poly_mulmid_rows_one_word): _nmod_poly_mulmid uses them for m
    output coefficients and len2 = 3 to 16 below 6 bits, 32 below 10 bits
    and 64 above, when m >= max(8, len2), in place of the rows (from len2 =
    6; below, only where the rows do not apply), of the classical code (one
    dot product per output, which has no one-word path), of KS and of
    fft_small; not for full products, which go to _nmod_poly_mul, nor for
    squares, where the classical code computes half of the products. With
    NEON, for m >= 16 and len2 up to 10 below 14 bits, 64 above, and before
    the rows. Measured (Nov 2026): within 3% of the best of the other
    algorithms from 8 outputs for len2 = 5 to 14 with gcc, 10 to 16 with
    clang, from about len2 outputs above (Cascade Lake, Ice Lake, Zen 4);
    on the M4 from 16 outputs, and KS faster for len2 >= 12 below 14 bits.
*/
static inline int
_nmod_poly_mulmid_want_one_word(slong len2, slong m, nmod_t mod)
{
    slong bits = NMOD_BITS(mod);

    if (_NMOD_POLY_NEON)
        return m >= 16 && len2 >= 3 && len2 <= (bits < 14 ? 10 : 64)
            && _nmod_poly_one_word(len2, mod);

    return m >= FLINT_MAX(8, len2) && len2 >= 3
        && len2 <= (bits < 6 ? 16 : bits < 10 ? 32 : 64)
        && _nmod_poly_one_word(len2, mod);
}

void _nmod_poly_mulmid_rows_one_word(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

/* the shortest row of _nmod_poly_mulmid_rows: the row length is a concave
   function of the index of the coefficient of poly2, so it is that of the
   first or the last coefficient used */
static inline slong
_nmod_poly_mulmid_rows_min_len(slong len1, slong len2, slong nlo, slong nhi)
{
    slong j0 = FLINT_MAX(0, nlo - len1 + 1);
    slong j1 = FLINT_MIN(len2, nhi) - 1;

    return FLINT_MIN(FLINT_MIN(j0 + len1, nhi) - FLINT_MAX(j0, nlo),
                     FLINT_MIN(j1 + len1, nhi) - FLINT_MAX(j1, nlo));
}

/*
    Classical products with AVX512-IFMA (_nmod_poly_mulmid_ifma, any window
    of the product, any n): for a shorter operand of
    _NMOD_POLY_IFMA_MIN_LEN(bits) to at most 511 coefficients, used in place
    of the rows, of the classical code, of KS and of fft_small. The
    function exists in every build (it runs the classical code without
    IFMA), the dispatchers only call it when NMOD_POLY_HAVE_IFMA.
    _nmod_poly_mulmid_ifma_crt is the same, but takes the variant by CRT
    for n > 2^52 whatever the lengths (for the tests).

    _NMOD_POLY_IFMA_MIN_LEN: below 8 coefficients the rows are faster
    (rows / kernel 0.82-0.86 at 7 on Ice Lake, 31 to 50 bits), except for
    52-bit moduli, where the rows run the slower scalar products of 51 to
    63 bits: the kernel from 4 (1.05-1.22x faster at 4, 1.5-1.9x at 6 and
    7; Zen 4 and Ice Lake, Nov 2026). Above 52 bits the kernel has larger
    fixed costs (the shorter operand split or reduced, the two-limb REDC):
    it starts from FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52 coefficients for
    53- to 63-bit moduli (Ice Lake 8, at least 0.99x the dispatchers
    without it; Zen 4 12: 0.93-0.95x at 10, 0.99-1.00x at 12), from 4 for
    64-bit moduli, which have no rows (1.5x or more at 4 on both), and is
    used for at least 16 output coefficients (nout), for odd n and powers
    of 2 (for even n, the scalar reduction of each output, NMOD_RED3, makes
    it 0.1-0.75x the dispatchers without it up to 16 coefficients).

    The largest len2 (above, fft_small) comes from flint-mparam.h:
    FLINT_NMOD_POLY_IFMA_MAX_LEN(bits) for full and middle products, the
    crossover with the dispatchers without the kernel on products of
    lengths L and L, 4L and L, 16L and L (the unbalanced ones cross first).
    It follows the number of 50-bit primes fft_small needs: one up to 20
    bits, two up to 45 bits (the kernel takes one IFMA per product up to 26
    bits, two above), three above (seven IFMA above 52 bits, or the CRT);
    Zen 4 crosses about twice as early as Ice Lake from 21 bits on. Modulo
    the 50-bit primes of fft_small, which need no CRT,
    FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME. Low products cross later: up to
    511 coefficients (at the primes of fft_small,
    FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME). Above 52 bits, the CRT
    variant from FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN (0: never): on Zen 4,
    where vpmullq is fast, it is 1.25-1.5x faster than seven IFMA from 64
    coefficients; on Ice Lake, 1.3-2x slower. The bound above 52 bits is
    256 on Ice Lake (448 for balanced products) and 128 on Zen 4 (with the
    CRT; 192 for balanced products).
*/
#if FLINT_BITS == 64 && defined(__AVX512F__) && defined(__AVX512DQ__) \
        && defined(__AVX512IFMA__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_POLY_HAVE_IFMA 1
#else
# define NMOD_POLY_HAVE_IFMA 0
#endif

#if !defined(FLINT_NMOD_POLY_IFMA_MAX_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME) \
    || !defined(FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME) \
    || !defined(FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN) \
    || !defined(FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52)
# error "flint-mparam.h must define FLINT_NMOD_POLY_IFMA_MAX_LEN, FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME, FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME, FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN and FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52"
#endif

#define _NMOD_POLY_IFMA_MIN_LEN(bits) \
    ((bits) < 52 ? 8 : (bits) == 52 || (bits) == 64 ? 4 \
     : FLINT_NMOD_POLY_IFMA_MIN_LEN_ABOVE_52)

/* whether n is one of the primes of the default fft_small context */
int _nmod_poly_is_fft_small_prime(nmod_t mod);

/* nout: the number of output coefficients; low: whether the product is a
   low product (mullow) */
static inline int
_nmod_poly_want_ifma(slong len2, slong nout, nmod_t mod, int low)
{
#if NMOD_POLY_HAVE_IFMA
    slong bits = NMOD_BITS(mod);

    if (len2 < _NMOD_POLY_IFMA_MIN_LEN(bits)
            || len2 > (low ? 511 : FLINT_NMOD_POLY_IFMA_MAX_LEN(bits)))
        return 0;

    if (mod.n > (UWORD(1) << 52))
        return nout >= 16 && ((mod.n & 1) || (mod.n & (mod.n - 1)) == 0);

    return len2 <= (low ? FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME
                        : FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME)
        || bits != 50 || !_nmod_poly_is_fft_small_prime(mod);
#else
    /* no kernel (and 2^52 does not fit in a 32-bit word) */
    return 0;
#endif
}

void _nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);
void _nmod_poly_mulmid_ifma_crt(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

void _nmod_poly_inv_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod);
void _nmod_poly_div_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr P, slong Plen, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod);
int nmod_poly_irreducible_binomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_trinomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_tetranomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_pentanomial(nmod_poly_t res, ulong n);
void _nmod_poly_divrem_q0_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, nn_srcptr B, slong lenA, ulong invL, nmod_t mod);

/* product with a precomputed transform, see geometric_progression.c */
struct mul_precomp_struct;
int _nmod_geometric_progression_mulmid_fft(nn_ptr z, slong zl, slong zh,
        nn_srcptr a, slong an, const struct mul_precomp_struct * F,
        slong s, slong bn, nmod_t mod);

#endif
