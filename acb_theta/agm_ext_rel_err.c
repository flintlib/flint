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
acb_theta_agm_ext_rel_err(arf_t err, const arf_t c2, const arf_t r,
                          slong nb_good, slong prec)
{
    fmpz_t exp;
    arb_t x, y, z;

    fmpz_init(exp);
    arb_init(x);
    arb_init(y);
    arb_init(z);
    
    arb_set_arf(x, r);
    arb_mul_2exp_si(x, x, 2);
    arb_sub_si(x, x, 1, prec);

    if (nb_good < 1) /* Need at least 1 good step */
    {
        arf_pos_inf(err);
    }
    else if (!arb_is_negative(x)) /* Convergence rate too slow */
    {
        arf_pos_inf(err);
    }
    else
    {
        fmpz_one(exp);
        fmpz_mul_2exp(exp, exp, FLINT_MAX(nb_good - 1, 0));
        arb_set_arf(x, r);
        arb_pow_fmpz(x, x, exp, prec); /* x = r^(2^(n-1)) */
        arb_mul_si(z, x, -2, prec);
        arb_add_si(z, z, 1, prec);
        arb_mul_arf(x, x, c2, prec); /* x = c2 r^(2^(n-1)) */
        arb_neg(y, x);
        arb_add_si(y, y, 1, prec);
        arb_mul(z, z, y, prec); /* z = (1-c2 r^(2^(n-1)))(1 - 2r^(2^(n-1))) */
        arb_mul_2exp_si(x, x, 1);
        arb_div(z, x, z, prec);
        
        arb_expm1(z, z, prec);
        arb_get_ubound_arf(err, z, prec);
    }

    fmpz_clear(exp);
    arb_clear(x);
    arb_clear(y);
    arb_clear(z);
}
