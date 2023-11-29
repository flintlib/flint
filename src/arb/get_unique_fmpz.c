/*
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_get_unique_fmpz(fmpz_t z, const arb_t x)
{
    if (!arb_is_finite(x))
    {
        return 0;
    }
    else if (arb_is_exact(x))
    {
        /* x = b*2^e, e >= 0 */
        if (arf_is_int(arb_midref(x)))
        {
            /* arf_get_fmpz aborts on overflow */
            arf_get_fmpz(z, arb_midref(x), ARF_RND_DOWN);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    /* if the radius is >= 1, there are at least two integers */
    else if (mag_cmp_2exp_si(arb_radref(x), 0) >= 0)
    {
        return 0;
    }
    /* there are 0 or 1 integers if the radius is < 1 */
    else
    {
        fmpz_t a, b, exp;
        int res;

        /* if the midpoint is exactly an integer, it is what we want */
        if (arf_is_int(arb_midref(x)))
        {
            /* arf_get_fmpz aborts on overflow */
            arf_get_fmpz(z, arb_midref(x), ARF_RND_DOWN);
            return 1;
        }

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(exp);

        /* if the radius is tiny, it can't be an integer */
        arf_bot(a, arb_midref(x));

        if (fmpz_cmp(a, MAG_EXPREF(arb_radref(x))) > 0)
        {
            res = 0;
        }
        else
        {
            arb_get_interval_fmpz_2exp(a, b, exp, x);

            if (COEFF_IS_MPZ(*exp))
            {
                flint_throw(FLINT_ERROR, "arb_get_unique_fmpz: input too large\n");
            }

            if (*exp >= 0)
            {
                res = fmpz_equal(a, b);

                if (res)
                {
                    fmpz_mul_2exp(a, a, *exp);
                    fmpz_mul_2exp(b, b, *exp);
                }
            }
            else
            {
                fmpz_cdiv_q_2exp(a, a, -(*exp));
                fmpz_fdiv_q_2exp(b, b, -(*exp));
                res = fmpz_equal(a, b);
            }

            if (res)
                fmpz_set(z, a);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(exp);

        return res;
    }
}

