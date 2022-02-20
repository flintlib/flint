/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

void
fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong exp)
{
    slong g1 = *g;
    slong u1;
    ulong bits;

    u1 = FLINT_ABS(*g);
    bits = FLINT_BIT_COUNT(u1);

    if (exp * bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        if (bits == 0)
            exp = (exp != 0);

L1:     /* Fits in small fmpz, even if exp == 0 and g is big */
        u1 = 1;
        for (; exp != 0; exp--)
            u1 *= g1;
L2:     if (COEFF_IS_MPZ(*f))
            _fmpz_clear_mpz(*f);
        *f = u1;
    }
    else if (u1 == 1)
    {
        exp &= 1;
        goto L1;
    }
    else if (bits <= SMALL_FMPZ_BITCOUNT_MAX) /* g is small */
    {
        __mpz_struct mg = {1, 1, NULL};
        __mpz_struct * mf;

        mg._mp_d = (mp_limb_t *) &u1;
        if (g1 < 0)
            mg._mp_size = -1;

        if (!COEFF_IS_MPZ(*f))
        {
            mf = _fmpz_new_mpz();
            *f = PTR_TO_COEFF(mf);
        }
        else
        {
            mf = COEFF_TO_PTR(*f);
        }

        flint_mpz_pow_ui(mf, &mg, exp);
        if ((mf->_mp_size == 1 || mf->_mp_size == -1) && mf->_mp_d[0] <= COEFF_MAX)
        {
            g1 = *f;
            *f = (mf->_mp_size > 0) ? mf->_mp_d[0] : -mf->_mp_d[0];
            _fmpz_clear_mpz(g1);
        }
    }
    else /* g is large and exp > 0 */
    {
        __mpz_struct * mf;
        if (!COEFF_IS_MPZ(*f))
        {
            mf = _fmpz_new_mpz();
            *f = PTR_TO_COEFF(mf);
        }
        else
            mf = COEFF_TO_PTR(*f);
        flint_mpz_pow_ui(mf, COEFF_TO_PTR(g1), exp);
    }
}
