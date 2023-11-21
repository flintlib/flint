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
arb_zeta_ui_vec_even(arb_ptr x, ulong start, slong num, slong prec)
{
    slong i;

    for (i = 0; i < num; i++)
        arb_zeta_ui(x + i, start + 2 * i, prec);
}

