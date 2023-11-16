/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_set_fmpz_2exp_fmpz(mag_t z, const fmpz_t man, const fmpz_t exp)
{
    if (fmpz_is_zero(man))
    {
        mag_zero(z);
    }
    else
    {
        mp_limb_t m;
        slong cexp;

        m = fmpz_abs_ubound_ui_2exp(&cexp, man, MAG_BITS);
        MAG_MAN(z) = m;
        _fmpz_add_fast(MAG_EXPREF(z), exp, cexp + MAG_BITS);
    }
}

void
mag_set_fmpz_2exp_fmpz_lower(mag_t z, const fmpz_t man, const fmpz_t exp)
{
    if (fmpz_is_zero(man))
    {
        mag_zero(z);
    }
    else
    {
        mp_limb_t m;
        slong cexp;

        m = fmpz_abs_lbound_ui_2exp(&cexp, man, MAG_BITS);
        MAG_MAN(z) = m;
        _fmpz_add_fast(MAG_EXPREF(z), exp, cexp + MAG_BITS);
    }
}

