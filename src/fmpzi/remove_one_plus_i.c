/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

slong
fmpzi_remove_one_plus_i(fmpzi_t res, const fmpzi_t x)
{
    slong s, t;
    int odd;

    if (fmpzi_is_zero(x))
    {
        fmpzi_zero(res);
        return 0;
    }

    if (fmpz_is_zero(fmpzi_imagref(x)))
    {
        s = fmpz_val2(fmpzi_realref(x));
        odd = 0;
    }
    else if (fmpz_is_zero(fmpzi_realref(x)))
    {
        s = fmpz_val2(fmpzi_imagref(x));
        odd = 0;
    }
    else
    {
        s = fmpz_val2(fmpzi_realref(x));
        t = fmpz_val2(fmpzi_imagref(x));

        if (s == t)
        {
            odd = 1;
        }
        else
        {
            s = FLINT_MIN(s, t);
            odd = 0;
        }
    }

    if (s == 0)
    {
        fmpzi_mul_i_pow_si(res, x, -s);
    }
    else
    {
        fmpz_tdiv_q_2exp(fmpzi_realref(res), fmpzi_realref(x), s);
        fmpz_tdiv_q_2exp(fmpzi_imagref(res), fmpzi_imagref(x), s);
        fmpzi_mul_i_pow_si(res, res, -s);
    }

    if (odd)
    {
        /* (a+bi) / (1+i) = ((a+b)/2) + ((b-a)/2)i */
        fmpz_t t;
        fmpz_init(t);
        fmpz_add(t, fmpzi_realref(res), fmpzi_imagref(res));
        fmpz_sub(fmpzi_imagref(res), fmpzi_imagref(res), fmpzi_realref(res));
        fmpz_tdiv_q_2exp(fmpzi_realref(res), t, 1);
        fmpz_tdiv_q_2exp(fmpzi_imagref(res), fmpzi_imagref(res), 1);
        fmpz_clear(t);
    }

    return 2 * s + odd;
}
