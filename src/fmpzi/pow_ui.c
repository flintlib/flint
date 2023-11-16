/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

/* todo: maybe handle (+-1 +- i)^n as a special case */

void
fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp)
{
    fmpzi_struct * R, * S, * T;
    fmpzi_t tmp;
    ulong bit;

    /* todo: cubing formula? */
    if (exp <= 2)
    {
        if (exp == 0)
            fmpzi_one(res);
        else if (exp == 1)
            fmpzi_set(res, x);
        else
            fmpzi_sqr(res, x);
        return;
    }

    if (fmpz_is_zero(fmpzi_imagref(x)))
    {
        fmpz_pow_ui(fmpzi_realref(res), fmpzi_realref(x), exp);
        fmpz_zero(fmpzi_imagref(res));
        return;
    }

    if (fmpz_is_zero(fmpzi_realref(x)))
    {
        fmpz_pow_ui(fmpzi_realref(res), fmpzi_imagref(x), exp);
        fmpz_zero(fmpzi_imagref(res));
        if (exp % 4 >= 2)
            fmpz_neg(fmpzi_realref(res), fmpzi_realref(res));
        if (exp % 2 == 1)
            fmpz_swap(fmpzi_realref(res), fmpzi_imagref(res));
        return;
    }

    if (res == x)
    {
        fmpzi_t tmp;
        fmpzi_init(tmp);
        fmpzi_pow_ui(tmp, x, exp);
        fmpzi_swap(tmp, res);
        fmpzi_clear(tmp);
        return;
    }

    fmpzi_init(tmp);

    R = res;
    S = tmp;

    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);

    fmpzi_sqr(R, x);

    if (bit & exp)
    {
        fmpzi_mul(S, R, x);
        T = R;
        R = S;
        S = T;
    }

    while (bit >>= 1)
    {
        fmpzi_sqr(S, R);

        if (bit & exp)
        {
            fmpzi_mul(R, S, x);
        }
        else
        {
            T = R;
            R = S;
            S = T;
        }
    }

    if (R != res)
        fmpzi_swap(R, res);

    fmpzi_clear(tmp);
}
