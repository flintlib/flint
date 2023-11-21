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
arb_zeta_ui_vec(arb_ptr x, ulong start, slong num, slong prec)
{
    slong i, num_odd, num_even, start_odd;
    arb_ptr tmp;

    num_odd = num / 2 + (start & num & 1);
    num_even = num - num_odd;

    start_odd = start % 2;

    arb_zeta_ui_vec_even(x, start + start_odd, num_even, prec);
    arb_zeta_ui_vec_odd(x + num_even, start + !start_odd, num_odd, prec);

    /* interleave */
    tmp = flint_malloc(sizeof(arb_struct) * num);
    for (i = 0; i < num_even; i++) tmp[i] = x[i];
    for (i = 0; i < num_odd; i++) tmp[num_even + i] = x[num_even + i];
    for (i = 0; i < num_even; i++) x[start_odd + 2  * i] = tmp[i];
    for (i = 0; i < num_odd; i++) x[!start_odd + 2  * i] = tmp[num_even + i];
    flint_free(tmp);
}

