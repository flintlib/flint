/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz.h"
#include "ulong_extras.h"

slong
fmpz_clog(const fmpz_t n, const fmpz_t b)
{
    fmpz_t t;
    int sign;
    slong r;

    if (fmpz_is_one(n))
        return 0;

    if (!COEFF_IS_MPZ(*b))
        return fmpz_clog_ui(n, *b);

    if (fmpz_cmp(n, b) <= 0)
        return 1;

    r = fmpz_dlog(n) / fmpz_dlog(b);

    fmpz_init(t);
    fmpz_pow_ui(t, b, r);
    sign = fmpz_cmp(t, n);

    /* Adjust down */
    if (sign > 0)
    {
        while (sign > 0)
        {
            fmpz_divexact(t, t, b);
            sign = fmpz_cmp(t, n);
            r--;
        }
        r += (sign != 0);
    }
    /* Adjust up */
    else if (sign < 0)
    {
        while (sign < 0)
        {
            fmpz_mul(t, t, b);
            sign = fmpz_cmp(t, n);
            r++;
        }
    }

    fmpz_clear(t);
    return r;
}

slong
fmpz_clog_ui(const fmpz_t n, ulong b)
{
    fmpz_t t;
    int sign;
    slong r;

    if (fmpz_is_one(n))
        return 0;

    if (b == 2)
    {
        fmpz_init(t);
        fmpz_sub_ui(t, n, 1);
        r = fmpz_bits(t);
        fmpz_clear(t);
        return r;
    }

    if (!COEFF_IS_MPZ(*n))
        return n_clog(*n, b);

    if (fmpz_cmp_ui(n, b) <= 0)
        return 1;

    r = fmpz_dlog(n) / log(b);

    fmpz_init(t);
    fmpz_set_ui(t, b);
    fmpz_pow_ui(t, t, r);
    sign = fmpz_cmp(t, n);

    /* Adjust down */
    if (sign > 0)
    {
        while (sign > 0)
        {
            fmpz_divexact_ui(t, t, b);
            sign = fmpz_cmp(t, n);
            r--;
        }
        r += (sign != 0);
    }
    /* Adjust up */
    else if (sign < 0)
    {
        while (sign < 0)
        {
            fmpz_mul_ui(t, t, b);
            sign = fmpz_cmp(t, n);
            r++;
        }
    }

    fmpz_clear(t);
    return r;
}
