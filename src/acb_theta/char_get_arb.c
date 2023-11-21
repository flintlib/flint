/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_theta.h"

void
acb_theta_char_get_arb(arb_ptr v, ulong a, slong g)
{
    slong k;

    for (k = g - 1; k >= 0; k--)
    {
        arb_set_si(&v[k], a & 1);
        a = a >> 1;
    }
    _arb_vec_scalar_mul_2exp_si(v, v, g, -1);
}
