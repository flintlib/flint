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
arb_zeta_ui_vec_odd(arb_ptr x, ulong start, slong num, slong prec)
{
    slong i, num_borwein;
    ulong cutoff;

    cutoff = 40 + 0.3 * prec;

    if (cutoff > start)
    {
        num_borwein = 1 + (cutoff - start) / 2;
        num_borwein = FLINT_MIN(num_borwein, num);
    }
    else
        num_borwein = 0;

    arb_zeta_ui_vec_borwein(x, start, num_borwein, 2, prec);
    for (i = num_borwein; i < num; i++)
        arb_zeta_ui(x + i, start + 2 * i, prec);
}

