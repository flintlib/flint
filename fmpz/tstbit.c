/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

int fmpz_tstbit(const fmpz_t f, ulong i)
{
    if (!COEFF_IS_MPZ(*f))
    {
        if (i < FLINT_BITS)
        {
            return ((WORD(1) << i) & *f) != 0;
        }
        else  /* i >= FLINT_BITS */
        {
            return *f < 0;
        }
    }
    else
    {
        return mpz_tstbit(COEFF_TO_PTR(*f), i);
    }
}

