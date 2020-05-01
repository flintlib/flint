/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

/* todo: might be worth it to detect p/q * root of unity */
/* todo: manually zero imaginary part? */
void
ca_qqbar_abs(ca_qqbar_t res, const ca_qqbar_t x)
{
    if (ca_qqbar_is_real(x))
    {
        if (ca_qqbar_real_sgn(x) >= 0)
            ca_qqbar_set(res, x);
        else
            ca_qqbar_neg(res, x);
    }
    else
    {
        ca_qqbar_t t;
        ca_qqbar_init(t);

        if (ca_qqbar_real_sgn(x) == 0)
        {
            ca_qqbar_i(t);
            ca_qqbar_mul(res, x, t);
            if (ca_qqbar_real_sgn(res) < 0)
                ca_qqbar_neg(res, res);
        }
        else
        {
            ca_qqbar_conj(t, x);
            ca_qqbar_mul(t, x, t);
            ca_qqbar_sqrt(res, t);
        }

        ca_qqbar_clear(t);
    }
}
