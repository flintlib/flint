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
acb_theta_ql_nb_steps(const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = n_clog(FLINT_MAX(1, prec), 2);
    
    if (g == 1)
    {
        return FLINT_MAX(0, n-8);
    }
    else
    {
        return FLINT_MAX(0, n-4);
    }
}
