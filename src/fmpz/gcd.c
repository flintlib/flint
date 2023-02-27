/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        ulong u1;

        if (c1 == 0)
        {
            fmpz_abs(f, h);
            return;
        }

        u1 = FLINT_ABS(c1);
        if (!COEFF_IS_MPZ(c2))  /* and h is also small */
        {
            ulong u2;

            if (c2 == 0)
            {
                fmpz_abs(f, g);
                return;
            }

            u2 = FLINT_ABS(c2);
            fmpz_set_ui(f, mpn_gcd_1((mp_srcptr) &u2, (mp_size_t) 1, u1));
        }
        else                    /* but h is large */
        {
            __mpz_struct * mpzc2 = COEFF_TO_PTR(c2);
            mp_size_t size = mpzc2->_mp_size;
            /* The sign is stored in the size of an mpz, and gcd_1 only takes
             * positive integers. */
            fmpz_set_ui(f, mpn_gcd_1(mpzc2->_mp_d, FLINT_ABS(size), u1));
        }
    }
    else                        /* g is large */
    {
        if (!COEFF_IS_MPZ(c2))  /* but h is small */
        {
            ulong u2;
            __mpz_struct * mpzc1;
            mp_size_t size;

            if (c2 == 0)
            {
                fmpz_abs(f, g);
                return;
            }

            u2 = FLINT_ABS(c2);
            mpzc1 = COEFF_TO_PTR(c1);
            size = mpzc1->_mp_size;
            fmpz_set_ui(f, mpn_gcd_1(mpzc1->_mp_d, FLINT_ABS(size), u2));
        }
        else
        {
            /* TODO: Change to mpn_gcd in order to save some calculations that
             * have already been already made. */
            __mpz_struct * mf = _fmpz_promote(f);
            mpz_gcd(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);
        }
    }
}
