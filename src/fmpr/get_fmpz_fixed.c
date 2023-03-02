/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

int
fmpr_get_fmpz_fixed_fmpz(fmpz_t y, const fmpr_t x, const fmpz_t exp)
{
    slong shift;

    if (fmpr_is_zero(x))
    {
        fmpz_zero(y);
        return 0;
    }

    shift = _fmpz_sub_small(fmpr_expref(x), exp);

    if (shift >= 0)
    {
        fmpz_mul_2exp(y, fmpr_manref(x), shift);
        return 0;
    }
    else
    {
        fmpz_tdiv_q_2exp(y, fmpr_manref(x), -shift);
        return 1;
    }
}

int
fmpr_get_fmpz_fixed_si(fmpz_t y, const fmpr_t x, slong exp)
{
    int result;
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, exp);
    result = fmpr_get_fmpz_fixed_fmpz(y, x, t);
    fmpz_clear(t);
    return result;
}
