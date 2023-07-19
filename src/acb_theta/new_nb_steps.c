/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_ql_new_nb_steps(const arb_mat_t cho, slong d, slong prec)
{
    slong lowprec = ACB_THETA_ELD_DEFAULT_PREC;
    arb_t x;
    slong res;

    arb_init(x);
    
    arb_set(x, arb_mat_entry(cho, d, d));
    arb_sqr(x, x, lowprec);
    arb_mul_2exp_si(x, x, -prec);
    arb_log(x, x, lowprec);
    /* This is dangerous. */
    res = arf_get_si(arb_midref(x), ARF_RND_NEAR) - 4;

    arb_clear(x);
    return res;
}
