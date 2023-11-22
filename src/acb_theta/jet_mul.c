/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static int
acb_theta_jet_le(const slong * tup1, const slong * tup2, slong g)
{
    slong k;

    for (k = 0; k < g; k++)
    {
        if (tup1[k] > tup2[k])
        {
            return 0;
        }
    }
    return 1;
}

void
acb_theta_jet_mul(acb_ptr res, acb_srcptr v1, acb_srcptr v2, slong ord, slong g, slong prec)
{
    slong nb = acb_theta_jet_nb(ord, g);
    acb_ptr aux;
    slong * tups;
    slong * diff;
    slong j, k, l;

    aux = _acb_vec_init(nb);
    tups = flint_malloc(nb * g * sizeof(slong));
    diff = flint_malloc(g * sizeof(slong));

    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < nb; k++)
        {
            if (!acb_theta_jet_le(tups + k * g, tups + j * g, g))
            {
                continue;
            }
            for (l = 0; l < g; l++)
            {
                diff[l] = tups[j * g + l] - tups[k * g + l];
            }
            acb_addmul(&aux[j], &v1[k], &v2[acb_theta_jet_index(diff, g)], prec);
        }
    }

    _acb_vec_set(res, aux, nb);

    _acb_vec_clear(aux, nb);
    flint_free(tups);
    flint_free(diff);
}
