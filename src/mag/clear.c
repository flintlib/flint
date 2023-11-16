/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_clear(mag_t x)
{
    /* fmpz_clear(), but avoids a redundant store */
    if (COEFF_IS_MPZ(MAG_EXP(x)))
        _fmpz_clear_mpz(MAG_EXP(x));
}

