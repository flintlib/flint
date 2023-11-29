/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

ulong
z_gcdinv(ulong * inv, slong a, ulong b)
{
    ulong g, ua = FLINT_ABS(a);

    if (ua >= b)
        ua %= b;

    g = n_gcdinv(inv, ua, b);

    if (a < WORD(0))
        *inv = n_submod(UWORD(0), *inv, b);

    return g;
}

int
fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;
    int val;

    if (fmpz_is_zero(h))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_invmod). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            ulong inv, gcd;
            if (c2 < WORD(0))
                c2 = -c2;
            if (c2 == WORD(1))
            {
                fmpz_zero(f);
                return 1;       /* special case not handled by n_invmod */
            }

            gcd = z_gcdinv(&inv, c1, c2);

            return (gcd == UWORD(1) ? fmpz_set_si(f, inv), 1 : 0);
        }
        else                    /* h is large and g is small */
        {
            __mpz_struct temp;  /* put g into a temporary mpz_t */
            __mpz_struct * mf;

            if (c1 < WORD(0))
            {
                c1 = -c1;
                temp._mp_d = (mp_limb_t *) & c1;
                temp._mp_size = -1;
            }
            else if (c1 == WORD(0))
                temp._mp_size = 0;
            else
            {
                temp._mp_d = (mp_limb_t *) & c1;
                temp._mp_size = 1;
            }

            mf = _fmpz_promote(f);
            val = mpz_invert(mf, &temp, COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* inverse mod h may result in small value */

            return val;
        }
    }
    else                        /* g is large */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is small */
        {
            ulong gcd, inv, r;
            if (c2 < WORD(0))
                c2 = -c2;
            if (c2 == WORD(1))
            {
                fmpz_zero(f);
                return 1;       /* special case not handled by z_gcd_invert */
            }
            /* reduce g mod h first */

            r = flint_mpz_fdiv_ui(COEFF_TO_PTR(c1), c2);

            gcd = z_gcdinv(&inv, r, c2);

            return (gcd == UWORD(1) ? fmpz_set_si(f, inv), 1 : 0);
        }
        else                    /* both are large */
        {
            __mpz_struct * mf = _fmpz_promote(f);
            val = mpz_invert(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* reduction mod h may result in small value */

            return val;
        }
    }
}
