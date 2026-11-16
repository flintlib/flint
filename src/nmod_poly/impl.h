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
