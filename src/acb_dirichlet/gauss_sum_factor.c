/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_gauss_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{

    slong k;
    acb_t tmp;

    /* early check */
    for (k = (G->neven == 2); k < G->num; k++)
    {
        /* if e > 1 and not primitive, 0 */
        if (chi->log[k] % G->P[k].p == 0 && G->P[k].e > 1)
        {
            acb_zero(res);
            return;
        }
    }

    /* factor */
    acb_one(res);
    acb_init(tmp);

    for (k = (G->neven == 2); k < G->num; k++)
    {
        ulong pe = G->P[k].pe.n;
        dirichlet_group_t Gp;
        dirichlet_char_t chip;

        dirichlet_subgroup_init(Gp, G, pe);
        dirichlet_char_init(chip, Gp);

        chip->n = chi->n % pe;

        if (k == 1 && G->neven == 2)
        {
            chip->log[0] = chi->log[0];
            chip->log[1] = chi->log[1];
        }
        else
            chip->log[0] = chi->log[k];

        /* chi_pe(a, q/pe) * G_pe(a) */
        acb_dirichlet_gauss_sum(tmp, Gp, chip, prec);
        acb_mul(res, res, tmp, prec);

        acb_dirichlet_chi(tmp, Gp, chip, (G->q / pe) % pe, prec);
        acb_mul(res, res, tmp, prec);

        dirichlet_char_clear(chip);
        dirichlet_group_clear(Gp);
    }

    if (G->q_even == 2)
        acb_neg(res, res);

    acb_clear(tmp);
}
