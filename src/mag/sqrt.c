/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

#ifdef __GNUC__
# define frexp __builtin_frexp
# define ldexp __builtin_ldexp
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

void
mag_sqrt(mag_t y, const mag_t x)
{
    if (mag_is_special(x))
    {
        mag_set(y, x);
    }
    else
    {
        double t;
        slong e;

        t = MAG_MAN(x) * ldexp(1.0, -MAG_BITS);

        if (MAG_IS_LAGOM(x))
        {
            e = MAG_EXP(x);

            if (e % 2 != 0)
            {
                e = (e - 1) >> 1;
                t *= 2.0;
            }
            else
            {
                e >>= 1;
            }

            t = sqrt(t) * (1 + 1e-13);
            _fmpz_demote(MAG_EXPREF(y));
            MAG_SET_D_2EXP(MAG_MAN(y), MAG_EXP(y), t, e);
        }
        else
        {
            if (fmpz_is_odd(MAG_EXPREF(x)))
                t *= 2.0;
            fmpz_fdiv_q_2exp(MAG_EXPREF(y), MAG_EXPREF(x), 1);
            t = sqrt(t) * (1 + 1e-13);
            mag_set_d_2exp_fmpz(y, t, MAG_EXPREF(y));
        }
    }
}

void
mag_sqrt_lower(mag_t y, const mag_t x)
{
    if (mag_is_special(x))
    {
        mag_set(y, x);
    }
    else
    {
        double t;
        slong e;

        t = MAG_MAN(x) * ldexp(1.0, -MAG_BITS);

        if (MAG_IS_LAGOM(x))
        {
            e = MAG_EXP(x);

            if (e % 2 != 0)
            {
                e = (e - 1) >> 1;
                t *= 2.0;
            }
            else
            {
                e >>= 1;
            }

            t = sqrt(t) * (1 - 1e-13);
            _fmpz_demote(MAG_EXPREF(y));
            MAG_SET_D_2EXP_LOWER(MAG_MAN(y), MAG_EXP(y), t, e);
        }
        else
        {
            if (fmpz_is_odd(MAG_EXPREF(x)))
                t *= 2.0;
            fmpz_fdiv_q_2exp(MAG_EXPREF(y), MAG_EXPREF(x), 1);
            t = sqrt(t) * (1 - 1e-13);
            mag_set_d_2exp_fmpz_lower(y, t, MAG_EXPREF(y));
        }
    }
}

