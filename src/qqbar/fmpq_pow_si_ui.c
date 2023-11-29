/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "qqbar.h"

void
qqbar_fmpq_pow_si_ui(qqbar_t res, const fmpq_t x, slong a, ulong b)
{
    if (b == 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (a == 0)
    {
        qqbar_one(res);
    }
    else if (a == 1)
    {
        qqbar_fmpq_root_ui(res, x, b);
    }
    else if (fmpq_sgn(x) >= 0)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpq_pow_si(t, x, a);
        qqbar_fmpq_root_ui(res, t, b);
        fmpq_clear(t);
    }
    else
    {
        qqbar_fmpq_root_ui(res, x, b);
        if (a > 0)
        {
            qqbar_pow_ui(res, res, a);
        }
        else
        {
            qqbar_inv(res, res);
            qqbar_pow_ui(res, res, -a);
        }
    }
}
