/*
    Copyright (C) 2012, 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_mul(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpz xv, yv;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        if (fmpr_is_zero(x) && fmpr_is_finite(y))
        {
            fmpr_zero(z);
        }
        else if (fmpr_is_zero(y) && fmpr_is_finite(x))
        {
            fmpr_zero(z);
        }
        else if ((fmpr_is_inf(x) && (fmpr_is_inf(y) || !fmpr_is_special(y))) ||
            (fmpr_is_inf(y) && !fmpr_is_special(x)))
        {
            if (fmpr_sgn(x) == fmpr_sgn(y))
                fmpr_pos_inf(z);
            else
                fmpr_neg_inf(z);
        }
        else
        {
            fmpr_nan(z);
        }

        return FMPR_RESULT_EXACT;
    }

    xv = *fmpr_manref(x);
    yv = *fmpr_manref(y);

    if (!COEFF_IS_MPZ(xv) && !COEFF_IS_MPZ(yv))
    {
        return _fmpr_mul_1x1(z, FLINT_ABS(xv), fmpr_expref(x),
            FLINT_ABS(yv), fmpr_expref(y), (xv ^ yv) < 0, prec, rnd);
    }
    else
    {
        slong xn, yn;
        int xsign, ysign;
        mp_limb_t xtmp, ytmp;
        mp_ptr xptr, yptr;

        FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, xv)
        FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, yv)

        if (xn >= yn)
            return _fmpr_mul_mpn(z, xptr, xn, fmpr_expref(x),
                yptr, yn, fmpr_expref(y), xsign ^ ysign, prec, rnd);
        else
            return _fmpr_mul_mpn(z, yptr, yn, fmpr_expref(y),
                xptr, xn, fmpr_expref(x), ysign ^ xsign, prec, rnd);
    }
}

