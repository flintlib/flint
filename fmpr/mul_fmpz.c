/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_mul_fmpz(fmpr_t z, const fmpr_t x, const fmpz_t y, slong prec, fmpr_rnd_t rnd)
{
    fmpz xv, yv;
    fmpz yexp;

    if (fmpr_is_special(x) || fmpz_is_zero(y))
    {
        if (fmpr_is_zero(x))
        {
            fmpr_zero(z);
        }
        else if (fmpz_is_zero(y) && fmpr_is_finite(x))
        {
            fmpr_zero(z);
        }
        else if (fmpr_is_inf(x) && !fmpz_is_zero(y))
        {
            if (fmpr_sgn(x) == fmpz_sgn(y))
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
    yv = *y;

    if (!COEFF_IS_MPZ(xv) && !COEFF_IS_MPZ(yv))
    {
        mp_limb_t ytmp;
        unsigned int bc;
        ytmp = FLINT_ABS(yv);
        count_trailing_zeros(bc, ytmp);
        ytmp >>= bc;
        yexp = bc;

        return _fmpr_mul_1x1(z, FLINT_ABS(xv), fmpr_expref(x),
            ytmp, &yexp, (xv ^ yv) < 0, prec, rnd);
    }
    else
    {
        slong xn, yn;
        int xsign, ysign;
        mp_limb_t xtmp, ytmp;
        mp_ptr xptr, yptr;
        yexp = 0;

        FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, xv)
        FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, yv)

        if (xn >= yn)
            return _fmpr_mul_mpn(z, xptr, xn, fmpr_expref(x),
                yptr, yn, &yexp, xsign ^ ysign, prec, rnd);
        else
            return _fmpr_mul_mpn(z, yptr, yn, &yexp,
                xptr, xn, fmpr_expref(x), ysign ^ xsign, prec, rnd);
    }
}

slong
fmpr_mul_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)
{
    fmpz xv;
    fmpz yexp;
    slong xn;
    int xsign, ysign;
    mp_limb_t xtmp, ytmp;
    mp_ptr xptr;

    if (fmpr_is_special(x) || (y == 0))
    {
        if (fmpr_is_zero(x))
        {
            fmpr_zero(z);
        }
        else if ((y == 0) && fmpr_is_finite(x))
        {
            fmpr_zero(z);
        }
        else if (fmpr_is_inf(x) && (y != 0))
        {
            if (fmpr_sgn(x) == ((y > 0) - (y < 0)))
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
    ytmp = FLINT_ABS(y);
    ysign = y < 0;
    yexp = 0;

    if (!COEFF_IS_MPZ(xv))
    {
        unsigned int bc;
        count_trailing_zeros(bc, ytmp);
        ytmp >>= bc;
        yexp = bc;

        return _fmpr_mul_1x1(z, FLINT_ABS(xv), fmpr_expref(x),
            ytmp, &yexp, (xv < 0) ^ ysign, prec, rnd);
    }
    else
    {
        FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, xv)

        return _fmpr_mul_mpn(z, xptr, xn, fmpr_expref(x),
            &ytmp, 1, &yexp, xsign ^ ysign, prec, rnd);
    }
}


slong
fmpr_mul_ui(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)
{
    fmpz xv;
    fmpz yexp;
    slong xn;
    int xsign;
    mp_limb_t xtmp, ytmp;
    mp_ptr xptr;

    if (fmpr_is_special(x) || (y == 0))
    {
        if (fmpr_is_zero(x))
        {
            fmpr_zero(z);
        }
        else if ((y == 0) && fmpr_is_finite(x))
        {
            fmpr_zero(z);
        }
        else if (fmpr_is_inf(x) && (y != 0))
        {
            if (fmpr_sgn(x) == (y != 0))
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
    ytmp = y;
    yexp = 0;

    if (!COEFF_IS_MPZ(xv))
    {
        unsigned int bc;
        count_trailing_zeros(bc, ytmp);
        ytmp >>= bc;
        yexp = bc;

        return _fmpr_mul_1x1(z, FLINT_ABS(xv), fmpr_expref(x),
            ytmp, &yexp, xv < 0, prec, rnd);
    }
    else
    {
        FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, xv)

        return _fmpr_mul_mpn(z, xptr, xn, fmpr_expref(x),
            &ytmp, 1, &yexp, xsign, prec, rnd);
    }
}
