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
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    __mpz_struct * mf;
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c1))
    {
        if (!COEFF_IS_MPZ(c2))
        {
            ulong th, tl;
            smul_ppmm(th, tl, c1, c2);
            fmpz_set_signed_uiui(f, th, tl);
            return;
        }
        else if (c1 != 0)
        {
            mf = _fmpz_promote(f);
            flint_mpz_mul_si(mf, COEFF_TO_PTR(c2), c1);
            return;
        }
    }
    
    if (!COEFF_IS_MPZ(*f))
    {
        if (c1 == 0 || c2 == 0)
        {
            *f = 0;
            return;
        }
        mf = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mf);
    }
    else
    {
        if (c1 == 0 || c2 == 0)
        {
            _fmpz_clear_mpz(*f);
            *f = 0;
            return;
        }
        
        mf = COEFF_TO_PTR(*f);
    }

    if (!COEFF_IS_MPZ(c2))
        flint_mpz_mul_si(mf, COEFF_TO_PTR(c1), c2);
    else
        mpz_mul(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
}
