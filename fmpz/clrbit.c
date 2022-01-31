/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_clrbit(fmpz_t f, ulong i)
{
    if (!COEFF_IS_MPZ(*f))
    {
        if (i < SMALL_FMPZ_BITCOUNT_MAX)
        {
            *f &= ~(WORD(1) << i);
        }
        /* i >= FLINT_BITS  --> nop */
    }
    else
    {
        __mpz_struct *ptr = COEFF_TO_PTR(*f);

        mpz_clrbit(ptr, i);
        _fmpz_demote_val(f);
    }
}

