/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x)
{
    if (!arb_is_finite(x))
    {
        flint_throw(FLINT_ERROR, "expected finite input in %s\n", __func__);
    }
    else if (arb_is_exact(x))
    {
        arf_get_fmpz_2exp(a, exp, arb_midref(x));
        fmpz_set(b, a);
    }
    else if (arf_is_zero(arb_midref(x)))
    {
        arf_t t;
        arf_init_set_mag_shallow(t, arb_radref(x));
        arf_get_fmpz_2exp(b, exp, t);
        fmpz_neg(a, b);
    }
    else
    {
        arf_t rad;
        fmpz_t tmp;
        slong shift;
        flint_bitcnt_t aval, bval;

        fmpz_init(tmp);

        arf_get_fmpz_2exp(a, exp, arb_midref(x));
        arf_init_set_mag_shallow(rad, arb_radref(x));
        arf_get_fmpz_2exp(b, tmp, rad);

        shift = _fmpz_sub_small(exp, tmp);

        if (FLINT_ABS(shift) >= WORD_MAX / 2)
            flint_throw(FLINT_ERROR, "too large shift in %s\n", __func__);

        if (shift >= 0)
        {
            fmpz_mul_2exp(a, a, shift);
            fmpz_set(exp, tmp);
        }
        else
            fmpz_mul_2exp(b, b, -shift);

        fmpz_sub(tmp, a, b);
        fmpz_add(b, a, b);
        fmpz_swap(tmp, a);

        if (fmpz_is_zero(a))
        {
            aval = fmpz_val2(b);
        }
        else if (fmpz_is_zero(b))
        {
            aval = fmpz_val2(a);
        }
        else
        {
            aval = fmpz_val2(a);
            bval = fmpz_val2(b);
            aval = FLINT_MIN(aval, bval);
        }

        if (aval > 0)
        {
            fmpz_add_ui(exp, exp, aval);
            fmpz_tdiv_q_2exp(a, a, aval);
            fmpz_tdiv_q_2exp(b, b, aval);
        }

        fmpz_clear(tmp);
    }
}

