/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

/* defined in exp.c */
void _mag_exp_d(mag_t res, double x, int roundup);
void mag_exp_huge(mag_t res, const mag_t x);
void mag_exp_huge_lower(mag_t res, const mag_t x);

void
mag_expinv(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
    }
    else if (mag_is_inf(x))
    {
        mag_zero(res);
    }
    else if (mag_cmp_2exp_si(x, 24) >= 0)
    {
        mag_exp_huge_lower(res, x);
        mag_inv(res, res);
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        mag_one(res);
    }
    else
    {
        slong e = MAG_EXP(x);

        if (e <= -MAG_BITS)
            mag_one(res);
        else
            _mag_exp_d(res, -ldexp(MAG_MAN(x), e - MAG_BITS), 1);
    }
}

void
mag_expinv_lower(mag_t res, const mag_t x)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
    }
    else if (mag_is_inf(x))
    {
        mag_zero(res);
    }
    else if (mag_cmp_2exp_si(x, 24) >= 0)
    {
        mag_exp_huge(res, x);
        mag_inv_lower(res, res);
    }
    else if (COEFF_IS_MPZ(MAG_EXP(x)))
    {
        /* 1 - eps */
        MAG_MAN(res) = (1 << MAG_BITS) - 1;
        fmpz_zero(MAG_EXPREF(res));
    }
    else
    {
        slong e = MAG_EXP(x);

        if (e < -MAG_BITS)
        {
            /* 1 - eps */
            MAG_MAN(res) = (1 << MAG_BITS) - 1;
            fmpz_zero(MAG_EXPREF(res));
        }
        else
            _mag_exp_d(res, -ldexp(MAG_MAN(x), e - MAG_BITS), 0);
    }
}

