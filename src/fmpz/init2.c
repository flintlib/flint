/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_init2(fmpz_t f, ulong limbs)
{
    if (limbs)
    {
        mpz_ptr mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        FLINT_MPZ_REALLOC(mf, limbs);
    }
    else
    {
        (*f) = WORD(0);
    }
}
