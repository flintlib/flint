/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

static void
_fmpr_mul_special(fmpr_t z, const fmpr_t x, const fmpr_t y)
{
    if (fmpr_is_zero(x))
    {
        if (!fmpr_is_special(y) || fmpr_is_zero(y))
            fmpr_zero(z);
        else
            fmpr_nan(z);
        return;
    }

    if (fmpr_is_zero(y))
    {
        if (!fmpr_is_special(x))
            fmpr_zero(z);
        else
            fmpr_nan(z);
        return;
    }

    if ((fmpr_is_inf(x) && (fmpr_is_inf(y) || !fmpr_is_special(y))) ||
        (fmpr_is_inf(y) && !fmpr_is_special(x)))
    {
        if (fmpr_sgn(x) == fmpr_sgn(y))
            fmpr_pos_inf(z);
        else
            fmpr_neg_inf(z);
        return;
    }

    fmpr_nan(z);
}

slong
fmpr_mul_naive(fmpr_t z, const fmpr_t x, const fmpr_t y, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        _fmpr_mul_special(z, x, y);
        return FMPR_RESULT_EXACT;
    }

    fmpz_mul(fmpr_manref(z), fmpr_manref(x), fmpr_manref(y));
    fmpz_add(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
    return _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, rnd);
}
