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

int
fmpz_sgn(const fmpz_t f)
{
    fmpz d = *f;

    if (d == 0)
        return 0;
    if (!COEFF_IS_MPZ(d))       /* c1 is small */
        return (d > WORD(0) ? 1 : -1);
    else
        return mpz_sgn(COEFF_TO_PTR(d));
}
