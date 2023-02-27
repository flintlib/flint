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
fmpz_set(fmpz_t f, const fmpz_t g)
{
    if (f == g)
        return;                 /* aliased inputs */

    if (!COEFF_IS_MPZ(*g))      /* g is small */
    {
        _fmpz_demote(f);
        *f = *g;
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);
        mpz_set(mf, COEFF_TO_PTR(*g));
    }
}
