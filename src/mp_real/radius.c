/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Predicates, magnitudes and radius updates. */

int
mp_real_is_zero(const mp_real_t x)
{
    return x->size == 0 && x->err == 0;
}

void
_mp_real_add_error_ulps(mp_real_t x, double e)
{
    _mp_real_apply_bnd(x, _mp_real_bnd_add(_mp_real_bnd_of(x),
        _mp_real_bnd_of_fberr(_mp_real_dbnd(e, x->exp - x->size))));
}

void
_mp_real_add_error_ulps_at(mp_real_t x, double v, slong anc)
{
    _mp_real_apply_bnd(x, _mp_real_bnd_add(_mp_real_bnd_of(x), _mp_real_bnd_of_fberr(_mp_real_dbnd(v, anc))));
}

void
mp_real_add_rel_error_2exp_si(mp_real_t x, slong e2)
{
    /* |value| < B^exp; add 2^e2 B^exp =
       2^(e2 mod FLINT_BITS) B^(exp + e2/FLINT_BITS), flooring the
       limb division so the bit remainder is >= 0 */
    slong q = e2 >> MP_REAL_LGB;
    int r = (int) (e2 - (q << MP_REAL_LGB));
    _mp_real_apply_bnd(x, _mp_real_bnd_add(_mp_real_bnd_of(x), _mp_real_bnd(0, UWORD(1) << r, x->exp + q)));
}

/* the relative radius of x (nonzero mantissa) is below 2^e: through
   the top limb, err < 2^bits(err); -WORD_MAX/2 when exact */
slong
mp_real_rel_radius_lt_2exp_si(const mp_real_t x)
{
    if (x->err == 0)
        return -WORD_MAX / 2;
    /* a zero midpoint with a radius: no relative accuracy at all */
    if (x->size == 0)
        return WORD_MAX / 2;
    return FLINT_BITS * (1 - x->size) + (slong) FLINT_BIT_COUNT(x->err)
        - (FLINT_BIT_COUNT(x->d[x->size - 1]) - 1) + 1;
}

/* x += [-2^e, 2^e] */
void
mp_real_add_error_2exp_si(mp_real_t x, slong e)
{
    slong q = e >> MP_REAL_LGB;
    _mp_real_add_error_ulps_at(x, (double) (UWORD(1) << (e - q * FLINT_BITS)), q);
}

/* e with |x| < 2^e for the ball (value plus radius), -WORD_MAX/2 for
   an exact zero */
slong
mp_real_abs_bound_lt_2exp_si(const mp_real_t x)
{
    slong e;

    if (x->size == 0)
        e = -WORD_MAX / 2;
    else
        e = FLINT_BITS * (x->exp - 1)
            + FLINT_BIT_COUNT(x->d[x->size - 1]);

    if (x->err != 0)
        e = FLINT_MAX(e, FLINT_BITS * (x->exp - x->size)
            + (slong) FLINT_BIT_COUNT(x->err)) + 1;
    return e;
}
