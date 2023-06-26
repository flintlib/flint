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
acb_theta_agm_nb_good_steps(const arf_t c, const arf_t r, slong prec)
{
    arb_t x;
    arb_t temp;
    arf_t b;
    fmpz_t exp;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb;

    arb_init(x);
    arb_init(temp);
    arf_init(b);
    fmpz_init(exp);

    /* Solve for c * r^(2^k) * (1+cr)/(1-r) <= 2^(-prec) */
    arb_one(x);
    arb_mul_2exp_si(x, x, -prec);
    arb_div_arf(x, x, c, lowprec);

    arb_set_arf(temp, r);
    arb_mul_arf(temp, temp, c, lowprec);
    arb_add_si(temp, temp, 1, lowprec);
    arb_div(x, x, temp, lowprec);

    arb_set_arf(temp, r);
    arb_sub_si(temp, temp, 1, lowprec);
    arb_neg(temp, temp);
    arb_mul(x, x, temp, lowprec);

    arb_log(x, x, lowprec);
    arb_set_arf(temp, r);
    arb_log(temp, temp, lowprec);
    arb_div(x, x, temp, lowprec);
    arb_get_ubound_arf(b, x, lowprec);

    if (arf_is_finite(b))
    {        
        arf_frexp(b, exp, b);
        nb = FLINT_MAX(0, fmpz_get_si(exp));
    }
    else
    {
        nb = -1; /* Error code */
    }

    arb_clear(x);
    arb_clear(temp);
    arf_clear(b);
    fmpz_clear(exp);
    return nb;
}
