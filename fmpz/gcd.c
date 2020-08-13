/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

ulong
z_gcd(slong a, slong b)
{
    ulong ua = FLINT_ABS(a);
    ulong ub = FLINT_ABS(b);

    return n_gcd(ua, ub);
}

void
fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(g))
    {
        fmpz_abs(f, h);
        return;
    }

    if (fmpz_is_zero(h))
    {
        fmpz_abs(f, g);
        return;
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz_set_si(f, z_gcd(c1, c2));
        }
        else                    /* h is large, but g is small */
        {
            fmpz c2d = fmpz_fdiv_ui(h, FLINT_ABS(c1));
            fmpz_set_si(f, z_gcd(c1, c2d));
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(c2))  /* h is small, but g is large */
        {
            fmpz c1d = fmpz_fdiv_ui(g, FLINT_ABS(c2));
            fmpz_set_si(f, z_gcd(c2, c1d));
        }
        else                    /* g and h are both large */
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* aliasing fine as g, h already large */

            mpz_gcd(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* gcd may be small */
        }
    }
}
