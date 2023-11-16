/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_pow_fmpz(mag_t z, const mag_t x, const fmpz_t e)
{
    if (fmpz_sgn(e) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, e);
        mag_inv(z, x);
        mag_pow_fmpz(z, z, t);
        fmpz_clear(t);
    }
    else if (!COEFF_IS_MPZ(*e))
    {
        mag_pow_ui(z, x, fmpz_get_ui(e));
    }
    else
    {
        mag_t y;
        mp_srcptr elimbs;
        slong i, bits;

        mag_init_set(y, x);
        bits = fmpz_bits(e);
        elimbs = COEFF_TO_PTR(*e)->_mp_d;

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul(y, y, y);

            if ((elimbs[i / FLINT_BITS] >> (i % FLINT_BITS)) & 1)
                mag_mul(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

void
mag_pow_fmpz_lower(mag_t z, const mag_t x, const fmpz_t e)
{
    if (fmpz_sgn(e) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, e);
        mag_inv_lower(z, x);
        mag_pow_fmpz_lower(z, z, t);
        fmpz_clear(t);
    }
    else if (!COEFF_IS_MPZ(*e))
    {
        mag_pow_ui_lower(z, x, fmpz_get_ui(e));
    }
    else
    {
        mag_t y;
        mp_srcptr elimbs;
        slong i, bits;

        mag_init_set(y, x);
        bits = fmpz_bits(e);
        elimbs = COEFF_TO_PTR(*e)->_mp_d;

        for (i = bits - 2; i >= 0; i--)
        {
            mag_mul_lower(y, y, y);

            if ((elimbs[i / FLINT_BITS] >> (i % FLINT_BITS)) & 1)
                mag_mul_lower(y, y, x);
        }

        mag_swap(z, y);
        mag_clear(y);
    }
}

