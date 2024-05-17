/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <gmp.h>
#include "fmpz.h"

void
fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    slong c1 = *g;
    ulong c1abs, c1bits;

    if (c1 == 0)
    {
        fmpz_zero(f);
        return;
    }

    c1abs = FLINT_ABS(c1);
    c1bits = FLINT_BIT_COUNT(c1abs);

    if (c1bits + exp <= SMALL_FMPZ_BITCOUNT_MAX)  /* Result fits inside a small fmpz */
    {
        if (COEFF_IS_MPZ(*f))
            _fmpz_clear_mpz(*f);

        *f = c1 << exp;
    }
    else if (c1bits <= SMALL_FMPZ_BITCOUNT_MAX)   /* g is small */
    {
        ulong expred = exp % FLINT_BITS;
        int alloc = 1 + exp / FLINT_BITS + ((c1bits + expred) > FLINT_BITS);
        ulong * limbs;
        mpz_ptr mf;

        mf = _fmpz_promote(f);
        limbs = FLINT_MPZ_REALLOC(mf, alloc);
        mf->_mp_size = (c1 > 0) ? alloc : -alloc;
        memset(limbs, 0, sizeof(ulong) * alloc);

        if (c1bits + expred <= FLINT_BITS)
        {
            limbs[alloc - 1] = c1abs << expred;
        }
        else
        {
            limbs[alloc - 1] = c1abs >> (FLINT_BITS - expred);
            limbs[alloc - 2] = c1abs << expred;
        }
    }
    else                                /* g is large */
    {
        mpz_ptr mg = COEFF_TO_PTR(c1);

        mpz_mul_2exp(_fmpz_promote(f), mg, exp);
    }
}
