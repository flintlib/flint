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
fmpr_set_fmpq(fmpr_t x, const fmpq_t y, slong prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_one(fmpq_denref(y)))
    {
        return fmpr_set_round_fmpz(x, fmpq_numref(y), prec, rnd);
    }
    else
    {
        slong res;

        fmpr_t t, u;
        fmpr_init(t);
        fmpr_init(u);

        fmpr_set_fmpz(t, fmpq_numref(y));
        fmpr_set_fmpz(u, fmpq_denref(y));

        res = fmpr_div(x, t, u, prec, rnd);

        fmpr_clear(t);
        fmpr_clear(u);

        return res;
    }
}
