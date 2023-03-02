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
fmpr_log(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_neg_inf(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else if (fmpr_sgn(x) < 0)
    {
        fmpr_nan(y);
        return FMPR_RESULT_EXACT;
    }
    else if (fmpr_is_one(x))
    {
        fmpr_zero(y);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong r;
        CALL_MPFR_FUNC(r, mpfr_log, y, x, prec, rnd);
        return r;
    }
}

slong
fmpr_log1p(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong r;
        CALL_MPFR_FUNC(r, mpfr_log1p, y, x, prec, rnd);
        return r;
    }
}
