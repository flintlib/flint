/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_vec.h"

void _arb_vec_set(arb_ptr res, arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_set(res + i, vec + i);
}

void _arb_vec_set_round(arb_ptr res, arb_srcptr vec, slong len, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_set_round(res + i, vec + i, prec);
}

void _arb_vec_zero(arb_ptr A, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_zero(A + i);
}

int _arb_vec_get_unique_fmpz_vec(fmpz * res,  arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!arb_get_unique_fmpz(res + i, vec + i))
            return 0;
    return 1;
}
