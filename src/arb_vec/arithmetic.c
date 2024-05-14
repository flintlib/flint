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

void _arb_vec_add(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_add(C + i, A + i, B + i, prec);
}

void _arb_vec_sub(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_sub(C + i, A + i, B + i, prec);
}

void _arb_vec_neg(arb_ptr B, arb_srcptr A, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_neg(B + i, A + i);
}
