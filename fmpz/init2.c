/*
    Copyright (C) 2009 William Hart

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
fmpz_init2(fmpz_t f, ulong limbs)
{
    if (limbs)
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        _mpz_realloc(mf, limbs);
    }
    else
    {
        (*f) = WORD(0);
    }
}
