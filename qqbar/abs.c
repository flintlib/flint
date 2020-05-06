/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

/* todo: might be worth it to detect p/q * root of unity */
/* todo: manually zero imaginary part? */
void
qqbar_abs(qqbar_t res, const qqbar_t x)
{
    if (qqbar_is_real(x))
    {
        if (qqbar_sgn_re(x) >= 0)
            qqbar_set(res, x);
        else
            qqbar_neg(res, x);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);

        if (qqbar_sgn_re(x) == 0)
        {
            qqbar_i(t);
            qqbar_mul(res, x, t);
            if (qqbar_sgn_re(res) < 0)
                qqbar_neg(res, res);
        }
        else
        {
            qqbar_conj(t, x);
            qqbar_mul(t, x, t);
            qqbar_sqrt(res, t);
        }

        qqbar_clear(t);
    }
}
