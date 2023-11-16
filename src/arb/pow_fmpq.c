/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, slong prec)
{
    if (fmpz_is_one(fmpq_denref(a)))
    {
        arb_pow_fmpz(y, x, fmpq_numref(a), prec);
    }
    else
    {
        int use_exp;
        slong k = *fmpq_denref(a);

        if (k == 2 || k == 4)
            use_exp = 0;
        else if (k > 1 && k < 50)
            use_exp = prec < (WORD(1) << ((k / 8) + 8));
        else
            use_exp = 1;

        if (use_exp)
        {
            arb_log(y, x, prec + 10);
            arb_mul_fmpz(y, y, fmpq_numref(a), prec + 10);
            arb_div_fmpz(y, y, fmpq_denref(a), prec + 10);
            arb_exp(y, y, prec);
        }
        else
        {
            arb_root_ui(y, x, k, prec);
            arb_pow_fmpz(y, y, fmpq_numref(a), prec);
        }
    }
}

