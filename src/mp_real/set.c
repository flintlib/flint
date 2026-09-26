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

/* Assignment and exact scaling. */

void
mp_real_set_ui(mp_real_t x, ulong c)
{
    if (c == 0)
    {
        mp_real_zero(x);
        return;
    }

    mp_real_fit_length(x, 1);
    x->d[0] = c;
    x->size = 1;
    x->negative = 0;
    x->exp = 1;
    x->err = 0;
}

void
mp_real_set_si(mp_real_t x, slong c)
{
    mp_real_set_ui(x, (c >= 0) ? (ulong) c : -(ulong) c);
    x->negative = (c < 0);
}

void
mp_real_set(mp_real_t res, const mp_real_t x)
{
    if (res == x)
        return;
    mp_real_fit_length(res, x->size);
    flint_mpn_copyi(res->d, x->d, x->size);
    res->size = x->size;
    res->negative = x->negative;
    res->exp = x->exp;
    res->err = x->err;
}

void
_mp_real_set_mpn_2exp(mp_real_t x, nn_srcptr p, slong len, slong ebits)
{
    slong b, q;

    while (len > 0 && p[len - 1] == 0)
        len--;
    if (len == 0)
    {
        /* keep the frame: a zero import anchors its ulp at
           B^ceil(ebits / FLINT_BITS) >= 2^ebits, so that a radius
           subsequently attached in ulps (a truncated-to-zero
           quantity known to |value| < k 2^ebits) lands at the
           intended scale instead of at B^0 */
        mp_real_zero(x);
        x->exp = ebits / FLINT_BITS + ((ebits % FLINT_BITS) > 0);
        return;
    }

    b = ebits % FLINT_BITS;
    if (b < 0)
        b += FLINT_BITS;
    q = (ebits - b) / FLINT_BITS;

    mp_real_fit_length(x, len + 1);
    if (b)
    {
        x->d[len] = mpn_lshift(x->d, p, len, (int) b);
        x->size = len + (x->d[len] != 0);
    }
    else
    {
        flint_mpn_copyi(x->d, p, len);
        x->size = len;
    }
    x->negative = 0;
    x->exp = x->size + q;
    x->err = 0;
    _mp_real_norm(x);
}

/* value *= 2^e (bit-level shift; everything else in mp_real is
   limb-aligned) */
void
mp_real_mul_2exp_si(mp_real_t res, const mp_real_t xin, slong e)
{
    mp_real_struct * x = res;

    if (res != xin)
        mp_real_set(res, xin);

    slong q = e >> MP_REAL_LGB;             /* floor limb division */
    int r = (int) (e - (q << MP_REAL_LGB)); /* 0 <= r < FLINT_BITS */

    if (x->size == 0)
    {
        x->exp += q + (r != 0);
        /* the radius scales exactly: err 2^r = err 2^(r - B) B at the
           anchor moved up by one limb when r != 0 */
        if (x->err != 0 && r != 0)
        {
            mp_real_bnd_t t = _mp_real_bnd(x->err >> (FLINT_BITS - r), x->err << r, x->exp - 1);
            _mp_real_zero_bnd(x, t);
        }
        return;
    }

    x->exp += q;
    if (r != 0)
    {
        ulong cy, err = x->err;

        mp_real_fit_length(x, x->size + 1);
        cy = mpn_lshift(x->d, x->d, x->size, r);
        if (cy)
        {
            x->d[x->size] = cy;
            x->size++;
            x->exp++;
        }
        _mp_real_norm(x);
        if (err != 0)
            _mp_real_apply_bnd(x, _mp_real_bnd(err >> (FLINT_BITS - r), err << r,
                x->exp - x->size));
        else
            x->err = 0;
    }
}
