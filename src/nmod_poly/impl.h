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
    a dot product on top of their common arithmetic. Not for powers of 2,
    whose dot products need no reduction (the rows lose from 52 bits).

    _NMOD_POLY_ROWS_RATIO: the ratio for full and low products, 4 with the
    AVX2 kernels above 32 bits (doubles, peeling to an aligned address,
    scalar tails), 2 with the 512-bit kernels (the only ones above 50 bits,
    masked tails) and up to 32 bits.

    _NMOD_POLY_ROWS_MAX_LEN2: the largest len2, 6 with the NEON kernels
    (on the M4 the rows are 1.15-2.0x faster for len2 = 3 to 6, and
    0.59-0.88x from 8, Nov 2026).
*/
#define _NMOD_POLY_ROWS_RATIO(mod) \
    ((NMOD_VEC_SIMD_SCALAR_MAX_BITS(1) > 50 || NMOD_BITS(mod) <= 32) ? 2 : 4)

#if defined(__ARM_NEON) && defined(__aarch64__)
# define _NMOD_POLY_ROWS_MAX_LEN2 6
#else
# define _NMOD_POLY_ROWS_MAX_LEN2 16
#endif

static inline int
_nmod_poly_want_rows(slong row_len, slong len2, slong ratio, nmod_t mod)
{
    return len2 >= 3 && len2 <= _NMOD_POLY_ROWS_MAX_LEN2
        && NMOD_BITS(mod) <= NMOD_VEC_SIMD_SCALAR_MAX_BITS(1)
        && row_len >= FLINT_MAX(NMOD_VEC_SIMD_MIN_LEN, ratio * len2)
        && (mod.n & (mod.n - 1)) != 0;
}

void _nmod_poly_mulmid_rows(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

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
    of the product): for moduli n <= 2^52 and a shorter operand of
    _NMOD_POLY_IFMA_MIN_LEN(bits) to at most 511 coefficients, used in place
    of the rows, of the classical code, of KS and of fft_small. The
    function exists in every build (it runs the classical code without
    IFMA), the dispatchers only call it when NMOD_POLY_HAVE_IFMA.

    _NMOD_POLY_IFMA_MIN_LEN: below 8 coefficients the rows are faster
    (rows / kernel 0.82-0.86 at 7 on Ice Lake, 31 to 50 bits), except for
    52-bit moduli, where the rows run the slower scalar products of 51 to
    63 bits: the kernel from 4 (1.05-1.22x faster at 4, 1.5-1.9x at 6 and
    7; Zen 4 and Ice Lake, Nov 2026).

    The largest len2 (above, fft_small) comes from flint-mparam.h:
    FLINT_NMOD_POLY_IFMA_MAX_LEN(bits) for full and middle products, the
    crossover with the dispatchers without the kernel on products of
    lengths L and L, 4L and L, 16L and L (the unbalanced ones cross first).
    It follows the number of 50-bit primes fft_small needs: one up to 20
    bits, two up to 45 bits (the kernel takes one IFMA per product up to 26
    bits, two above), three above; Zen 4 crosses about twice as early as
    Ice Lake from 21 bits on. Modulo the 50-bit primes of fft_small, which
    need no CRT, FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME. Low products cross
    later: up to 511 (at the primes of fft_small,
    FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME).
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
    || !defined(FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME)
# error "flint-mparam.h must define FLINT_NMOD_POLY_IFMA_MAX_LEN, FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME and FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME"
#endif

#define _NMOD_POLY_IFMA_MIN_LEN(bits) ((bits) == 52 ? 4 : 8)

/* whether n is one of the primes of the default fft_small context */
int _nmod_poly_is_fft_small_prime(nmod_t mod);

/* low: whether the product is a low product (mullow) */
static inline int
_nmod_poly_want_ifma(slong len2, nmod_t mod, int low)
{
    slong bits = NMOD_BITS(mod);

    if (!NMOD_POLY_HAVE_IFMA || mod.n > (UWORD(1) << 52)
            || len2 < _NMOD_POLY_IFMA_MIN_LEN(bits)
            || len2 > (low ? 511 : FLINT_NMOD_POLY_IFMA_MAX_LEN(bits)))
        return 0;

    return len2 <= (low ? FLINT_NMOD_POLY_IFMA_LOW_MAX_LEN_FFT_PRIME
                        : FLINT_NMOD_POLY_IFMA_MAX_LEN_FFT_PRIME)
        || bits != 50 || !_nmod_poly_is_fft_small_prime(mod);
}

void _nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod);

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
