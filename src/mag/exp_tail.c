/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_exp_tail(mag_t z, const mag_t x, ulong N)
{
    if (N == 0 || mag_is_inf(x))
    {
        mag_exp(z, x);
    }
    else if (mag_is_zero(x))
    {
        mag_zero(z);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_set_ui_2exp_si(t, N, -1);

        /* bound by geometric series when N >= 2*x  <=> N/2 >= x */
        if (mag_cmp(t, x) >= 0)
        {
            /* 2 c^N / N! */
            mag_pow_ui(t, x, N);
            mag_rfac_ui(z, N);
            mag_mul(z, z, t);
            mag_mul_2exp_si(z, z, 1);
        }
        else
        {
            mag_exp(z, x);
        }

        mag_clear(t);
    }
}
