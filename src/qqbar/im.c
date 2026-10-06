/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"
#include "impl.h"

void
qqbar_im(qqbar_t res, const qqbar_t x)
{
    if (qqbar_sgn_im(x) == 0)
    {
        qqbar_zero(res);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);

        if (qqbar_sgn_re(x) == 0)
        {
            qqbar_i(t);
            qqbar_mul(res, x, t);
            qqbar_neg(res, res);
        }
        else
        {
            /* im(x) = (x - conj(x)) / (2i) */
            qqbar_conj(t, x);

            if (qqbar_degree(x) >= 9 && _qqbar_binary_op_structured(res, x, t, 6))
            {
                /* conj(x) is in Q(x); we have computed im(x)^2 */
                qqbar_sqrt(res, res);
                if (qqbar_sgn_im(x) < 0)
                    qqbar_neg(res, res);
            }
            else
            {
                _qqbar_conjugate_pair_op(res, x, t, 5);
            }
        }

        arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));
        qqbar_clear(t);
    }
}
