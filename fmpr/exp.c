/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_exp(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_one(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else if (fmpr_is_neg_inf(x))
            fmpr_zero(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong r;
        CALL_MPFR_FUNC(r, mpfr_exp, y, x, prec, rnd);
        return r;
    }
}

slong
fmpr_expm1(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else if (fmpr_is_neg_inf(x))
            fmpr_set_si(y, -WORD(1));
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong r;
        CALL_MPFR_FUNC(r, mpfr_expm1, y, x, prec, rnd);
        return r;
    }
}
