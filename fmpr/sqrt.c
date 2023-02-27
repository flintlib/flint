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
fmpr_sqrt(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    slong r;

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

    if (fmpr_sgn(x) < 0)
    {
        fmpr_nan(y);
        return FMPR_RESULT_EXACT;
    }

    /* special case: 4^n */
    /* TODO: process all small exact square roots efficiently */
    if (fmpz_is_one(fmpr_manref(x)) && fmpz_is_even(fmpr_expref(x)))
    {
        r = fmpr_set_round(y, x, prec, rnd);
        fmpz_tdiv_q_2exp(fmpr_expref(y), fmpr_expref(y), 1);
        return r;
    }

    {
        fmpr_t t;
        fmpz_t e;

        fmpr_init(t);
        fmpz_init(e);

        fmpz_neg(e, fmpr_expref(x));
        if (fmpz_is_odd(e))
            fmpz_add_ui(e, e, 1);
        fmpr_mul_2exp_fmpz(t, x, e);

        CALL_MPFR_FUNC(r, mpfr_sqrt, y, t, prec, rnd);

        fmpz_neg(e, e);
        fmpz_tdiv_q_2exp(e, e, 1);
        fmpr_mul_2exp_fmpz(y, y, e);

        fmpr_clear(t);
        fmpz_clear(e);

        return r;
    }
}

