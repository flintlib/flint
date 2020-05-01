/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_im(ca_qqbar_t res, const ca_qqbar_t x)
{
    if (ca_qqbar_sgn_im(x) == 0)
    {
        ca_qqbar_zero(res);
    }
    else
    {
        ca_qqbar_t t;
        ca_qqbar_init(t);

        if (ca_qqbar_sgn_re(x) == 0)
        {
            ca_qqbar_i(t);
            ca_qqbar_mul(res, x, t);
            ca_qqbar_neg(res, res);
        }
        else
        {
            ca_qqbar_conj(t, x);
            ca_qqbar_sub(res, x, t);
            ca_qqbar_i(t);
            ca_qqbar_mul(res, res, t);
            ca_qqbar_neg(res, res);
            ca_qqbar_mul_2exp_si(res, res, -1);
        }

        arb_zero(acb_imagref(CA_QQBAR_ENCLOSURE(res)));
        ca_qqbar_clear(t);
    }
}
