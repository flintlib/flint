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
fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c2 = *g;

    if (x == 0)
    {
        fmpz_zero(f);
        return;
    }
    else if (!COEFF_IS_MPZ(c2)) /* c2 is small */
    {
        mp_limb_t prod[2];
        mp_limb_t uc2 = FLINT_ABS(c2);

        /* unsigned limb by limb multiply (assembly for most CPU's) */
        umul_ppmm(prod[1], prod[0], uc2, x);
        if (c2 < WORD(0))
            fmpz_neg_uiui(f, prod[1], prod[0]);
        else
            fmpz_set_uiui(f, prod[1], prod[0]);
    }
    else                        /* c2 is large */
    {
        /* Promote without val as if aliased both are large */
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        flint_mpz_mul_ui(mpz_ptr, COEFF_TO_PTR(c2), x);
    }
}
