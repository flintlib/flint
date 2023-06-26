/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
arb_randtest_pos(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    int pos = 0;
    while (!pos)
    {
        arb_randtest_precise(x, state, prec, mag_bits);
        pos = arb_is_positive(x);
    }
}
