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

void
acb_theta_jet_exp_pi_i(acb_ptr res, arb_srcptr a, slong ord, slong g, slong prec)
{
    slong nb = acb_theta_jet_nb(ord, g);
    slong * tups;
    acb_t c;
    arb_t t;
    fmpz_t den, m;
    slong k, l;

    tups = flint_malloc(nb * g * sizeof(slong));
    acb_init(c);
    arb_init(t);
    fmpz_init(den);
    fmpz_init(m);

    acb_one(&res[0]);
    acb_theta_jet_tuples(tups, ord, g);

    for (k = 1; k < nb; k++)
    {
        acb_one(&res[k]);
        fmpz_one(den);
        for (l = 0; l < g; l++)
        {
            arb_pow_ui(t, &a[l], tups[k * g + l], prec);
            acb_mul_arb(&res[k], &res[k], t, prec);

            fmpz_fac_ui(m, tups[k * g + l]);
            fmpz_mul(den, den, m);
        }

        acb_const_pi(c, prec);
        acb_mul_onei(c, c);
        acb_pow_ui(c, c, acb_theta_jet_total_order(tups + k * g, g), prec);
        acb_mul(&res[k], &res[k], c, prec);
        acb_div_fmpz(&res[k], &res[k], den, prec);
    }

    flint_free(tups);
    acb_clear(c);
    arb_clear(t);
    fmpz_clear(den);
    fmpz_clear(m);
}
