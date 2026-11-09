/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal definitions of acf. */

#ifndef ACF_IMPL_H
#define ACF_IMPL_H

#include "acf.h"
#include "../arf/elem.h"

/* the working precision beyond which the retries of the strict modes
   of the compositions give up and round their approximations as they
   stand (acf/elem.c) */
#define ACF_ELEM_RETRY_LIMIT_BITS 65536

/* the most bits of the spread of the input's exponents that the retries
   of the functions computed from acb balls add to their limit */
#define ACF_ELEM_SPREAD_MAX (WORD(1) << 20)

/* the directed modes of a value to be negated */
FLINT_FORCE_INLINE arf_rnd_t
_acf_rnd_flip(arf_rnd_t rnd)
{
    return (rnd == ARF_RND_FLOOR) ? ARF_RND_CEIL : (rnd == ARF_RND_CEIL) ? ARF_RND_FLOOR : rnd;
}

/* the cap of _acf_part_gap */
#define ACF_PART_GAP_MAX (WORD(1) << 30)

/* the exponent gap E(Rref) - E(R) >= 0 between a part R and the larger
   part Rref of the same complex number (0 when either is zero; acf/elem.c) */
slong _acf_part_gap(const arf_t R, const arf_t Rref);

/* the finite arf R of at most p bits, known within K ulps at p bits of
   R (the errors of the compositions) plus 2^ebu ulps at p bits of the
   larger part of the complex number (WORD_MIN for none), D the gap of
   _acf_part_gap, or exactly (exact = 1), rounded into res in the mode
   rnd by _arf_set_mpn_err: the inexact flag, or ARF_ELEM_UNCERTAIN when
   the window does not decide the rounding or the relaxed modes need a
   smaller error (the caller retries at a higher precision).  The
   relaxed modes measure the errors against an ulp of the larger part
   (the complex-relative accuracy of acf.h).  A zero R is an exact zero
   unless an absolute error is given (acf/elem.c). */
int _acf_elem_set(arf_t res, const arf_t R, slong D, slong p, ulong K, slong ebu, int exact, slong prec, arf_rnd_t rnd);

/* the exact square root of a nonzero x into res, or 0 when it is not
   dyadic (acf/elem.c) */
int _acf_sqrt_exact(acf_t res, const acf_t x);

/* Exact rounding decisions (acf/side.c): lists of exact terms, the sign
   of their sum for any exponents, the rounding boundary nearest to an
   approximation, and the rounding of a value from its side of that
   boundary. */
#define ACF_TERMS_MAX 32

typedef struct
{
    arf_struct t[ACF_TERMS_MAX];
    slong n;
}
_acf_terms_struct;

typedef _acf_terms_struct _acf_terms_t[1];

void _acf_terms_init(_acf_terms_t T);
void _acf_terms_clear(_acf_terms_t T);
void _acf_terms_add_monomial(_acf_terms_t T, slong c, arf_srcptr x1, arf_srcptr x2, arf_srcptr x3, arf_srcptr x4);
void _acf_terms_mul(_acf_terms_t R, const _acf_terms_t X, const _acf_terms_t Y, slong c);
void _acf_terms_append(_acf_terms_t R, const _acf_terms_t X);
int _acf_sign_sum(const _acf_terms_t T);
void _acf_round_boundary(arf_t beta, const arf_t y, slong prec, arf_rnd_t rnd);
int _acf_round_side(arf_t res, const arf_t beta, int side, slong prec, arf_rnd_t rnd);

/* the functions of _acf_near_side (acf/near.c): the side of a part of
   f(x) against a boundary, from the expansion at the nearest exact point */
enum
{
    ACF_NEAR_EXP, ACF_NEAR_EXP_PI_I, ACF_NEAR_SIN, ACF_NEAR_COS,
    ACF_NEAR_SIN_PI, ACF_NEAR_COS_PI, ACF_NEAR_LOG, ACF_NEAR_TAN,
    ACF_NEAR_COT, ACF_NEAR_SEC, ACF_NEAR_CSC, ACF_NEAR_TAN_PI,
    ACF_NEAR_COT_PI, ACF_NEAR_SEC_PI, ACF_NEAR_CSC_PI,
    ACF_NEAR_EXPM1, ACF_NEAR_LOG1P, ACF_NEAR_SINC, ACF_NEAR_SINC_PI,
    ACF_NEAR_ASIN, ACF_NEAR_ACOS, ACF_NEAR_ATAN, ACF_NEAR_ACOT,
    ACF_NEAR_ASEC, ACF_NEAR_ACSC, ACF_NEAR_ASINH, ACF_NEAR_ACOSH,
    ACF_NEAR_ATANH, ACF_NEAR_ACOTH, ACF_NEAR_ASECH, ACF_NEAR_ACSCH,
    ACF_NEAR_ASIN_PI, ACF_NEAR_ACOS_PI, ACF_NEAR_ATAN_PI, ACF_NEAR_NONE
};

int _acf_near_side(int * side, int fn, int part, const acf_t x, const arf_t beta, int bpi);

/* x^y near its exact cases (x = 2^K i^u (1 + h) with h tiny) and for a
   tiny y: 1 with the result and its flags, or 0 (acf/near.c) */
int _acf_pow_near(acf_t res, int * status, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd);

#endif
