/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_set_round_ui_2exp_fmpz(fmpr_t z,
        mp_limb_t lo, const fmpz_t exp, int negative,
        slong prec, fmpr_rnd_t rnd)
{
    slong lead, trail, bc, shift, shift2, ret;

    shift = 0;

    if ((lo & 1) == 0)
    {
        if (lo == 0)
        {
            fmpr_zero(z);
            return FMPR_RESULT_EXACT;
        }

        count_trailing_zeros(trail, lo);
        lo >>= trail;
        shift = trail;
    }

    count_leading_zeros(lead, lo);
    bc = FLINT_BITS - lead;

    ret = FMPR_RESULT_EXACT;
    if (bc > prec)
    {
        shift2 = bc - prec;
        lo = (lo >> shift2) + rounds_up(rnd, negative);
        count_trailing_zeros(trail, lo);
        lo >>= trail;
        shift += shift2;
        shift += trail;
        ret = trail;

        /* special case: if the mantissa overflowed to the next power of two,
           the error bound must be multiplied by two */
        ret -= (trail == prec);
    }

    if (!negative)
        fmpz_set_ui(fmpr_manref(z), lo);
    else
        fmpz_neg_ui(fmpr_manref(z), lo);

    fmpz_add_si_inline(fmpr_expref(z), exp, shift);
    return ret;
}

