/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz c1;

    if (exp == WORD(0))
    {
        fmpz_one(f);
        return;
    }

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        ulong u1 = FLINT_ABS(c1);
        ulong bits = FLINT_BIT_COUNT(u1);
        if (u1 <= UWORD(1))
        {
            fmpz_set_ui(f, u1);
        }
        else if (exp * bits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            fmpz_set_ui(f, n_pow(u1, exp));
        }
        else
        {
            __mpz_struct * mf = _fmpz_promote_val(f);

            flint_mpz_set_ui(mf, u1);
            flint_mpz_pow_ui(mf, mf, exp);
            _fmpz_demote_val(f);    /* may actually fit into a small after all */
        }

        if ((c1 < WORD(0)) && (exp & 1)) /* sign is -ve if exp odd and g -ve */
            fmpz_inplace_neg(f);
    }
    else
    {
        __mpz_struct * mf = _fmpz_promote_val(f);
        flint_mpz_pow_ui(mf, COEFF_TO_PTR(c1), exp);
        /* no need to demote as it can't get smaller */
    }
}

int fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e)
{
    int e_sgn = fmpz_sgn(e);

    if (e_sgn < 0)
    {
        flint_throw(FLINT_ERROR, "Negative exponent in fmpz_pow_fmpz");
    }
    else if (e_sgn == 0)
    {
        fmpz_one(a);
    }
    else if (fmpz_is_zero(b))
    {
        fmpz_zero(a);
    }
    else if (fmpz_is_pm1(b))
    {
        fmpz_set_si(a, fmpz_is_one(b) || fmpz_is_even(e) ? 1 : -1);
    }
    else
    {
        if (!fmpz_fits_si(e))
            return 0;

        fmpz_pow_ui(a, b, fmpz_get_si(e));
    }
    return 1;
}

void
fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e)
{
    if (e <= 1)
    {
        fmpz_set_ui(x, e == 0 ? 1 : b);
    }
    else if (e == 2)
    {
        mp_limb_t t[2];
        umul_ppmm(t[1], t[0], b, b);
        fmpz_set_uiui(x, t[1], t[0]);
    }
    else if (b <= 1)
    {
        fmpz_set_ui(x, b);
    }
    else
    {
        ulong bits = FLINT_BIT_COUNT(b);

        if (e * bits <= FLINT_BITS)
        {
            fmpz_set_ui(x, n_pow(b, e));
        }
        else
        {
            __mpz_struct * z = _fmpz_promote(x);
            flint_mpz_set_ui(z, b);
            flint_mpz_pow_ui(z, z, e);
            _fmpz_demote_val(x);
        }
    }
}
