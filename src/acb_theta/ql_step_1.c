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
acb_theta_ql_step_1(acb_ptr r, acb_srcptr th, acb_srcptr th0, acb_srcptr roots,
    arb_srcptr dist, arb_srcptr dist0, slong g, slong prec)
{
    slong n = 1 << g;
    
    flint_printf("(ql_step) input:");
    _acb_vec_printd(th, n, 10);
    flint_printf("\n");

    /* Be more careful with precisions */
    if (th == th0)
    {
        acb_theta_agm_sqr(r, th0, g, prec);        
    }
    else
    {
        acb_theta_agm_mul(r, th, th0, g, prec);
    }
    _acb_vec_scalar_mul_2exp_si(r, r, n, g);
    
    flint_printf("(ql_step) after duplication:\n");
    _acb_vec_printd(r, n, 10);
    flint_printf("\n");        
    flint_printf("(ql_step) square roots:\n");
    _acb_vec_printd(roots, n, 10);
    flint_printf("\n");
    
    acb_theta_agm_sqrt(r, r, roots, n, prec);
}
