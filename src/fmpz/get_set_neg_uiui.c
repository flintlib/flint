/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

void fmpz_get_uiui(ulong * hi, ulong * low, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        *low = *f;
        *hi  = 0;
    }
    else
    {
        mpz_ptr mpz = COEFF_TO_PTR(*f);
        *low = mpz->_mp_d[0];
        *hi  = mpz->_mp_size == 2 ? mpz->_mp_d[1] : 0;
    }
}

void fmpz_set_uiui(fmpz_t f, ulong hi, ulong lo)
{
    if (hi == 0)
    {
        fmpz_set_ui(f, lo);
    }
    else
    {
        mpz_ptr z = _fmpz_promote(f);
        if (z->_mp_alloc < 2)
            mpz_realloc2(z, 2 * FLINT_BITS);
        z->_mp_d[0] = lo;
        z->_mp_d[1] = hi;
        z->_mp_size = 2;
    }
}

void fmpz_neg_uiui(fmpz_t f, ulong hi, ulong lo)
{
    if (hi == 0)
    {
        fmpz_neg_ui(f, lo);
    }
    else
    {
        mpz_ptr z = _fmpz_promote(f);
        if (z->_mp_alloc < 2)
            mpz_realloc2(z, 2 * FLINT_BITS);
        z->_mp_d[0] = lo;
        z->_mp_d[1] = hi;
        z->_mp_size = -2;
    }
}
