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
acb_theta_agm_radius(arf_t rad, const arf_struct* mi, const arf_struct* Mi,
    const arf_t abs_dist, slong nb, slong prec)
{
    arb_t rho;
    arb_t next;
    arb_t t1, t2;
    slong k;

    arb_init(rho);
    arb_init(next);
    arb_init(t1);
    arb_init(t2);

    arb_set_arf(rho, abs_dist);
    for (k = 0; k < nb; k++)
    {
        /* Set rho to solution of sqrt((M+rho')/(m-rho')) * rho' <= rho */
        arb_set_arf(next, &mi[nb - 1 - k]);
        arb_mul_2exp_si(next, next, -1);
        arb_min(next, next, rho, prec);
        arb_add_arf(t1, next, &Mi[nb - 1 - k], prec);
        arb_sub_arf(t2, next, &mi[nb - 1 - k], prec);
        arb_neg(t2, t2);
        arb_div(t1, t1, t2, prec);
        arb_sqrt(t1, t1, prec);
        arb_div(next, rho, t1, prec);
        arb_set(rho, next);
    }

    arb_get_lbound_arf(rad, rho, prec);

    arb_clear(rho);
    arb_clear(next);
    arb_clear(t1);
    arb_clear(t2);
}
