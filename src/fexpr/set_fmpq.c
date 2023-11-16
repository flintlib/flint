/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_set_fmpq(fexpr_t res, const fmpq_t x)
{
    if (fmpz_is_one(fmpq_denref(x)))
    {
        fexpr_set_fmpz(res, fmpq_numref(x));
    }
    else
    {
        slong p, q;

        p = *fmpq_numref(x);
        q = *fmpq_denref(x);

        if (p >= FEXPR_COEFF_MIN && p <= FEXPR_COEFF_MAX &&
            q >= FEXPR_COEFF_MIN && q <= FEXPR_COEFF_MAX)
        {
            fexpr_fit_size(res, 4);
            res->data[0] = FEXPR_TYPE_CALL2 | (4 << FEXPR_TYPE_BITS);
            res->data[1] = FEXPR_SYMBOL_Div;
            res->data[2] = p << FEXPR_TYPE_BITS;
            res->data[3] = q << FEXPR_TYPE_BITS;
        }
        else   /* todo: copy directly without temporaries */
        {
            fexpr_t a, b;

            fexpr_init(a);
            fexpr_init(b);

            fexpr_set_fmpz(a, fmpq_numref(x));
            fexpr_set_fmpz(b, fmpq_denref(x));
            fexpr_div(res, a, b);

            fexpr_clear(a);
            fexpr_clear(b);
        }
    }
}
