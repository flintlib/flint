/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_mag_div_ui_ui(mag_t res, ulong a, ulong b)
{
    mag_set_ui(res, a);
    mag_div_ui(res, res, b);
}

void
acb_dirichlet_backlund_s_bound(mag_t res, const arb_t t)
{
    if (!arb_is_nonnegative(t))
    {
        mag_inf(res);
    }
    else
    {
        mag_t m;
        mag_init(m);
        arb_get_mag(m, t);
        if (mag_cmp_2exp_si(m, 8) < 0) /* 2^8 < 280 */
        {
            mag_one(res);
        }
        else if (mag_cmp_2exp_si(m, 22) < 0) /* 2^22 < 6.8*10^6 */
        {
            mag_set_ui(res, 2);
        }
        else if (mag_cmp_2exp_si(m, 29) < 0) /* 2^29 < 5.45*10^8 */
        {
            _mag_div_ui_ui(res, 231366, 100000);
        }
        else
        {
            /* |S(t)| <= 0.112*log(t) + 0.278*log(log(t)) + 2.51 */
            mag_t c, logm;
            mag_init(c);
            mag_init(logm);
            mag_log(logm, m);
            _mag_div_ui_ui(c, 278, 1000);
            mag_log(res, logm);
            mag_mul(res, res, c);
            _mag_div_ui_ui(c, 112, 1000);
            mag_addmul(res, c, logm);
            _mag_div_ui_ui(c, 251, 100);
            mag_add(res, res, c);
            mag_clear(c);
            mag_clear(logm);
        }
        mag_clear(m);
    }
}
