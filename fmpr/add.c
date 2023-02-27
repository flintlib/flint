/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
_fmpr_add_1x1(fmpr_t z,
        mp_limb_t x, int xsign, const fmpz_t xexp,
        mp_limb_t y, int ysign, const fmpz_t yexp,
        slong shift, slong prec, slong rnd)
{
    mp_limb_t hi, lo, t, u;
    int sign = ysign;

    t = x;
    u = y;

    (void) yexp; /* unused */

    lo = u << shift;
    hi = (shift == 0) ? 0 : (u >> (FLINT_BITS - shift));

    if (xsign == ysign)
    {
        add_ssaaaa(hi, lo, hi, lo, 0, t);
    }
    else
    {
        if (hi == 0)
        {
            if (lo >= t)
            {
                lo = lo - t;
            }
            else
            {
                lo = t - lo;
                sign = !sign;
            }
        }
        else
        {
            sub_ddmmss(hi, lo, hi, lo, 0, t);
        }
    }

    if (hi == 0)
        return fmpr_set_round_ui_2exp_fmpz(z, lo, xexp, sign, prec, rnd);
    else
        return fmpr_set_round_uiui_2exp_fmpz(z, hi, lo, xexp, sign, prec, rnd);
}

static slong
_fmpr_add_special(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_zero(x))
    {
        if (fmpr_is_zero(y))
        {
            fmpr_zero(z);
            return FMPR_RESULT_EXACT;
        }
        else
            return fmpr_set_round(z, y, prec, rnd);
    }
    else if (fmpr_is_zero(y))
    {
        return fmpr_set_round(z, x, prec, rnd);
    }
    else if (fmpr_is_nan(x) || fmpr_is_nan(y)
        || (fmpr_is_pos_inf(x) && fmpr_is_neg_inf(y))
        || (fmpr_is_neg_inf(x) && fmpr_is_pos_inf(y)))
    {
        fmpr_nan(z);
        return FMPR_RESULT_EXACT;
    }
    else if (fmpr_is_special(x))
    {
        fmpr_set(z, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        fmpr_set(z, y);
        return FMPR_RESULT_EXACT;
    }
}

slong
fmpr_add(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    slong shift, xn, yn;
    mp_limb_t xtmp, ytmp;
    mp_ptr xptr, yptr;
    fmpz xv, yv;
    const fmpz * xexp;
    const fmpz * yexp;
    int xsign, ysign;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        return _fmpr_add_special(z, x, y, prec, rnd);
    }

    shift = _fmpz_sub_small(fmpr_expref(y), fmpr_expref(x));

    if (shift >= 0)
    {
        xexp = fmpr_expref(x);
        yexp = fmpr_expref(y);
        xv = *fmpr_manref(x);
        yv = *fmpr_manref(y);
    }
    else
    {
        xexp = fmpr_expref(y);
        yexp = fmpr_expref(x);
        xv = *fmpr_manref(y);
        yv = *fmpr_manref(x);
        shift = -shift;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, xv)
    FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, yv)

    if ((xn == 1) && (yn == 1) && (shift < FLINT_BITS))
        return _fmpr_add_1x1(z, xptr[0], xsign, xexp, yptr[0], ysign, yexp, shift, prec, rnd);
    else
        return _fmpr_add_mpn(z, xptr, xn, xsign, xexp, yptr, yn, ysign, yexp, shift, prec, rnd);
}


