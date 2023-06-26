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
acb_theta_agm_conv_rate(arf_t c, arf_t r, const arf_t eps, slong prec)
{
    arb_t e;
    arb_t t;

    arb_init(t);
    arb_init(e);

    arb_set_arf(e, eps);
    arb_mul_2exp_si(e, e, 2); /* 4eps */
    arb_sub_si(e, e, 1, prec);

    if (!arb_is_negative(e))
    {
        arf_pos_inf(c);
        arf_pos_inf(r);
    }
    else
    {   
        /* Set t = eta(eps) */  
        arb_set_arf(e, eps);
        arb_add_si(e, e, 2, prec);
        arb_mul_arf(e, e, eps, prec); /* 2eps + eps^2 */
        
        arb_sub_si(t, e, 1, prec);
        arb_neg(t, t);
        arb_sqrt(t, t, prec);
        
        arb_mul_2exp_si(e, e, -1); /* eps + eps^2/2 */
        arb_add(t, t, e, prec);
        arb_neg(t, t);
        arb_add_si(t, t, 1, prec);
        
        arb_set_arf(e, eps);
        arb_sqr(e, e, prec); /* eps^2 */        
        arb_div(t, t, e, prec);
        
        arb_one(e);
        arb_mul_2exp_si(e, e, -1); /* 1/2 */
        arb_add(t, t, e, prec);
        
        arb_set_arf(e, eps);
        arb_sub_si(e, e, 1, prec);
        arb_neg(e, e); /* 1-eps */
        arb_div(t, t, e, prec);
        
        /* Set e to t*eps, t to 1/t */
        arb_mul_arf(e, t, eps, prec);
        arb_inv(t, t, prec);
        arb_get_ubound_arf(c, t, prec);
        arb_get_ubound_arf(r, e, prec);
    }

    arb_clear(t);
    arb_clear(e);
}
