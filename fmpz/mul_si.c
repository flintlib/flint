/*
    Copyright (C) 2020 William Hart
    Copyright (C) 2021 Albin Ahlbäck

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

void
fmpz_mul_si(fmpz_t f, const fmpz_t g, slong x)
{
    fmpz c2 = *g;

    if (!COEFF_IS_MPZ(c2)) /* c2 is small */
    {
        mp_limb_t th, tl;

        /* limb by limb multiply (assembly for most CPU's) */
        smul_ppmm(th, tl, c2, x);
        fmpz_set_signed_uiui(f, th, tl);
    }
    else                        /* c2 is large */
    {
        __mpz_struct * mf;
        if (!COEFF_IS_MPZ(*f))
        {
            if (x == 0)
            {
                *f = 0;
                return;
            }
            
            mf = _fmpz_new_mpz();
            *f = PTR_TO_COEFF(mf);
        }
        else
        {
            if (x == 0)
            {
                _fmpz_clear_mpz(*f);
                *f = 0;
                return;
            }

            mf = COEFF_TO_PTR(*f);
        }

        flint_mpz_mul_si(mf, COEFF_TO_PTR(c2), x);
    }
}
