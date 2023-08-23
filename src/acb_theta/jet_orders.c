/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_jet_orders(slong* orders, slong ord, slong g)
{
    slong nb_max, nb_rec;
    slong* rec;
    slong k, j, i, ind;

    if (g == 1)
    {
        orders[0] = ord;
        return;
    }

    nb_max = acb_theta_jet_nb(ord, g - 1);
    rec = flint_malloc(nb_max * (g - 1) * sizeof(slong));

    ind = 0;
    for (k = 0; k <= ord; k++)
    {
        nb_rec = acb_theta_jet_nb(k, g - 1);
        acb_theta_jet_orders(rec, k, g - 1);
        for (j = 0; j < nb_rec; j++)
        {
            orders[ind] = ord - k;
            for (i = 0; i < g - 1; i++)
            {
                orders[ind + 1 + i] = rec[j * (g - 1) + i];
            }
            ind += g;
        }
    }

    flint_free(rec);
}
