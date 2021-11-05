/*
    Copyright (C) 2020 William Hart

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
        mp_limb_t prod[2];
        mp_limb_t mc2 = c2;
        mp_limb_t mx = x;

        /* limb by limb multiply (assembly for most CPU's) */
        smul_ppmm(prod[1], prod[0], mc2, mx);
        fmpz_set_signed_uiui(f, prod[1], prod[0]);
    }
    else                        /* c2 is large */
    {
        if (x == 0)
            fmpz_zero(f);
        else
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* ok without val as if aliased both are large */
            flint_mpz_mul_si(mpz_ptr, COEFF_TO_PTR(c2), x);
        }
    }
}
