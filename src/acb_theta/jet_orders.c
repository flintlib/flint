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
acb_theta_jet_orders(slong* orders, slong ord, slong g)
{
    slong k, j, l, nb_rec, ind;
    slong* rec;

    if (g == 1)
    {
        for (k = 0; k <= ord; k++)
        {
            orders[k] = k;
        }
        return;
    }

    /* Generate orders in dimension g - 1 */
    nb_rec = acb_theta_jet_nb(ord, g - 1);
    rec = flint_malloc((g - 1) * nb_rec * sizeof(slong));
    acb_theta_jet_orders(rec, ord, g - 1);

    for (k = 0; k <= ord; k++)
    {
        /* Construct orders of total order k from rec */
        ind = acb_theta_jet_nb(k - 1, g);
        for (j = 0; j < acb_theta_jet_nb(k, g - 1); j++)
        {
            orders[(ind + j) * g] = k - acb_theta_jet_total_order(rec + j * (g - 1), g - 1);
            for (l = 0; l < g - 1; l++)
            {
                orders[(ind + j) * g + l + 1] = rec[j * (g - 1) + l];
            }
        }
    }

    flint_free(rec);
}
