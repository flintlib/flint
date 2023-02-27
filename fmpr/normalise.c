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
_fmpr_normalise_naive(fmpz_t man, fmpz_t exp, slong prec, fmpr_rnd_t rnd)
{
    /* TODO: this should perhaps raise an exception to avoid ambiguity */
    if (fmpz_is_zero(man))
    {
        fmpz_zero(exp);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong bc, val;

        bc = fmpz_bits(man);
        val = fmpz_val2(man);

        if (bc - val <= prec)
        {
            if (val != 0)
            {
                fmpz_tdiv_q_2exp(man, man, val);
                fmpz_add_ui(exp, exp, val);
            }

            return FMPR_RESULT_EXACT;
        }
        else
        {
            slong exp_shift = bc - prec;

            if (rnd == FMPR_RND_NEAR)
            {
                flint_abort();
            }
            else if (rnd == FMPR_RND_DOWN)
            {
                fmpz_tdiv_q_2exp(man, man, exp_shift);
            }
            else if (rnd == FMPR_RND_FLOOR)
            {
                fmpz_fdiv_q_2exp(man, man, exp_shift);
            }
            else if (rnd == FMPR_RND_CEIL)
            {
                fmpz_cdiv_q_2exp(man, man, exp_shift);
            }
            else
            {
                if (fmpz_sgn(man) > 0)
                    fmpz_cdiv_q_2exp(man, man, exp_shift);
                else
                    fmpz_fdiv_q_2exp(man, man, exp_shift);
            }

            val = fmpz_val2(man);
            exp_shift += val;

            if (val != 0)
                fmpz_tdiv_q_2exp(man, man, val);

            fmpz_add_ui(exp, exp, exp_shift);
            return val - (val == prec);
        }
    }
}
