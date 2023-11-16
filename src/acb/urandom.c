/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_urandom(acb_t z, flint_rand_t state, slong prec)
{
    arb_t abs;
    slong k;
    int done = 0;

    arb_init(abs);

    while (!done)
    {
        arb_urandom(acb_realref(z), state, prec);
        arb_urandom(acb_imagref(z), state, prec);
        /* Re-sample if z is outside the unit circle */
        acb_abs(abs, z, prec);
        arb_sub_si(abs, abs, 1, prec);
        done = arb_is_nonpositive(abs);
    }

    k = n_randint(state, 4);
    acb_mul_i_pow_si(z, z, k);

    arb_clear(abs);
}
