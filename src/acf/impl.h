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

#endif
