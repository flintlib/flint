/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "mag.h"

void
mag_set_d_2exp_fmpz(mag_t z, double c, const fmpz_t exp)
{
    if (c < 0.0)
        c = -c;

    if (c == 0.0)
    {
        mag_zero(z);
    }
    else if ((c != c) || c == D_INF)
    {
        mag_inf(z);
    }
    else
    {
        slong cexp = *exp;

        if (cexp >= MAG_MIN_LAGOM_EXP && cexp <= MAG_MAX_LAGOM_EXP)
        {
            _fmpz_demote(MAG_EXPREF(z));
            MAG_SET_D_2EXP(MAG_MAN(z), MAG_EXP(z), c, cexp);
        }
        else
        {
            MAG_SET_D_2EXP(MAG_MAN(z), cexp, c, 0);
            fmpz_add_si(MAG_EXPREF(z), exp, cexp);
        }
    }
}

void
mag_set_d_2exp_fmpz_lower(mag_t z, double c, const fmpz_t exp)
{
    if (c < 0.0)
        c = -c;

    if (c == 0.0 || (c != c))
    {
        mag_zero(z);
    }
    else if (c == D_INF)
    {
        mag_inf(z);
    }
    else
    {
        slong cexp = *exp;

        if (cexp >= MAG_MIN_LAGOM_EXP && cexp <= MAG_MAX_LAGOM_EXP)
        {
            _fmpz_demote(MAG_EXPREF(z));
            MAG_SET_D_2EXP_LOWER(MAG_MAN(z), MAG_EXP(z), c, cexp);
        }
        else
        {
            MAG_SET_D_2EXP_LOWER(MAG_MAN(z), cexp, c, 0);
            fmpz_add_si(MAG_EXPREF(z), exp, cexp);
        }
    }
}

