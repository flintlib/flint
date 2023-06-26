/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_deriv_get_dz(acb_t v, const acb_theta_deriv_t D, slong ab, slong ord, slong* indices)
{
    slong ord = 0;
    slong index = 0;
    slong g = acb_theta_deriv_dim(D);
    slong k;

    for (k = 0; k < g; k++)
    {
        ord += orders[k];
    }

    if (g == 1)
    {
        index = ;
    }
    else
    {
        index = (n_pow(g, ord) - 1) / (g - 1);
        for (k = 0; k < g; k++)
        {
            
        }
    }
    
    index = ab * (1 + binom(1, g) + binom(2, g+1) + ... + binom(acb_theta_deriv_ord(D), g + acb_theta_deriv_ord(D)-1)); /* this is binom(acb_theta_deriv_ord(D), g + acb_theta_deriv_ord(D)) */
    index += (1 + binom(1, g) + ... + binom(ord-1, g+ord-2)); /* this is binom, too */
    index += binom(ord, g-1+ord-1) + binom(ord-1, g-1+ord-2) + ... + binom(ord-orders[0]), bla)
    acb_set(v, acb_theta_deriv_val(D, k));
    
}

