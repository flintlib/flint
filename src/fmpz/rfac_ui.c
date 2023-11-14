/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

static inline ulong rfac(ulong x, ulong b)
{
    ulong i, c = x;

    for (i = 1; i < b; i++)
        c *= x + i;

    return c;
}

/* Assumes x positive, b > a. b must also be small enough to
avoid integer overflow, which is no problem if the result
is to fit in memory. */
void
_fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong a, ulong b)
{
    if (b - a == 1)
    {
        fmpz_add_ui(r, x, a);
    }
    else if ((*x <= COEFF_MAX) && (b - a < 60))
    {
        ulong step, bits, factors_per_limb;
        ulong y = *x;

        /* Bound size of largest factor */
        bits = FLINT_BIT_COUNT(y + a + b - 1);

        /* The result fits in a single limb */
        if ((b - a) * bits < FLINT_BITS)
            step = factors_per_limb = b - a;
        else
        {
            factors_per_limb = FLINT_BITS / bits;
            step = FLINT_MIN(b - a, factors_per_limb);
        }

        fmpz_set_ui(r, rfac(y + a, step));
        a += step;

        while (a < b)
        {
            step = FLINT_MIN(b - a, factors_per_limb);
            fmpz_mul_ui(r, r, rfac(y + a, step));
            a += step;
        }
    }
    else
    {
        fmpz_t t, u;
        ulong m = (a + b) / 2;

        fmpz_init(t);
        fmpz_init(u);

        _fmpz_rfac_ui(t, x, a, m);
        _fmpz_rfac_ui(u, x, m, b);
        fmpz_mul(r, t, u);

        fmpz_clear(t);
        fmpz_clear(u);
    }
}

void
fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong n)
{
    if (n == 0)
    {
        fmpz_one(r);
    }
    else if (n == 1)
    {
        fmpz_set(r, x);
    }
    else if (fmpz_is_zero(x))
    {
        fmpz_zero(r);
    }
    else if (fmpz_sgn(x) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_add_ui(t, x, n - 1);
        if (fmpz_sgn(t) >= 0)
        {
            fmpz_zero(r);
        }
        else
        {
            fmpz_neg(t, t);
            fmpz_rfac_ui(r, t, n);
            if (n % 2 == 1)
                fmpz_neg(r, r);
        }
        fmpz_clear(t);
    }
    else
    {
        _fmpz_rfac_ui(r, x, 0, n);
    }
}
