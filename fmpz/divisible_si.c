/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2022 Albin Ahlbäck

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

int
fmpz_divisible_si(const fmpz_t f, slong g)
{
    slong c = *f;

    if (c == WORD(0))
    {
        return 1;
    }

    if (!COEFF_IS_MPZ(c))
    {
        return !(c % g);
    }
    else
    {
        return flint_mpz_divisible_ui_p(COEFF_TO_PTR(c), FLINT_ABS(g));
    }
}

