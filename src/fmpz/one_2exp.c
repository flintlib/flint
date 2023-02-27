/*
    Copyright (C) 2021 Fredrik Johansson

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
fmpz_one_2exp(fmpz_t x, ulong e)
{
    if (e <= FLINT_BITS - 3)
    {
        fmpz_set_ui(x, UWORD(1) << e);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(x);
        flint_mpz_set_ui(z, 1);
        mpz_mul_2exp(z, z, e);
    }
}
