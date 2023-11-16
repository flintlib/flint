/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

/* todo: might be worth it to detect p/q * root of unity */
void
qqbar_abs2(qqbar_t res, const qqbar_t x)
{
    if (qqbar_is_real(x))
    {
        qqbar_sqr(res, x);
    }
    else if (qqbar_is_root_of_unity(NULL, NULL, x))
    {
        qqbar_one(res);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);

        if (qqbar_sgn_re(x) == 0)
        {
            qqbar_i(t);
            qqbar_mul(res, x, t);
            qqbar_sqr(res, res);
        }
        else
        {
            qqbar_conj(t, x);
            qqbar_mul(res, x, t);
        }

        qqbar_clear(t);
    }

    arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
}
