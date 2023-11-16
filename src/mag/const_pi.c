/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_const_pi(mag_t res)
{
    fmpz_clear(MAG_EXPREF(res));
    MAG_MAN(res) = 843314857;
    MAG_EXP(res) = 2;
}

void
mag_const_pi_lower(mag_t res)
{
    fmpz_clear(MAG_EXPREF(res));
    MAG_MAN(res) = 843314856;
    MAG_EXP(res) = 2;
}

