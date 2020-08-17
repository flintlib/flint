/*
    Copyright (C) 2011 Fredrik Johansson

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
fmpz_is_square(const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c <= 1)
            return (c >= 0);

        return n_is_square(c);
    }
    else
        return mpz_perfect_square_p(COEFF_TO_PTR(c));
}
