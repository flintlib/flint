/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        flint_printf("fmpz_mpoly_q_div_fmpz: division by zero\n");
        flint_abort();
    }
    else
    {
        if (fmpz_sgn(y) > 0)
        {
            fmpz_t one;
            *one = 1;
            _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                        fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                        one, y,
                        ctx);
        }
        else
        {
            fmpz_t t;
            fmpz_t one;
            *one = -1;
            fmpz_init(t);
            fmpz_neg(t, y);

            _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                        fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                        one, t,
                        ctx);

            fmpz_clear(t);
        }
    }

}

