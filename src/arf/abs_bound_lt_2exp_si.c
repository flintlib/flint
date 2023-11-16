/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

slong
arf_abs_bound_lt_2exp_si(const arf_t x)
{
    slong res;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            return -ARF_PREC_EXACT;
        else
            return ARF_PREC_EXACT;
    }

    if (!COEFF_IS_MPZ(ARF_EXP(x)))
        return ARF_EXP(x);

    if (fmpz_fits_si(ARF_EXPREF(x)))
        res = fmpz_get_si(ARF_EXPREF(x));
    else
        res = fmpz_sgn(ARF_EXPREF(x)) < 0 ? -ARF_PREC_EXACT : ARF_PREC_EXACT;

    if (res < -ARF_PREC_EXACT)
        res = -ARF_PREC_EXACT;
    if (res > ARF_PREC_EXACT)
        res = ARF_PREC_EXACT;

    return res;
}

