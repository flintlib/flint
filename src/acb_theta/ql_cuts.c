/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_ql_cuts(slong* cuts, const arb_mat_t cho, slong prec)
{
    slong k;
    slong g = arb_mat_nrows(cho);
    slong nb_cuts = 0;
    arb_t cmp;

    arb_init(cmp);

    for (k = 1; k < g; k++)
    {
        arb_mul_si(cmp, arb_mat_entry(cho, k - 1, k - 1),
            ACB_THETA_QL_CUT, prec);
        if (arb_lt(cmp, arb_mat_entry(cho, k, k)))
        {
            cuts[nb_cuts] = k;
            nb_cuts += 1;
        }
    }

    return nb_cuts;
}
