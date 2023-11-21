/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_jet_compose(acb_ptr res, acb_srcptr v, const acb_mat_t N,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(N);
    slong nb = acb_theta_jet_nb(ord, g);
    acb_ptr aux;
    acb_t x;
    fmpz_t m, p;
    slong * tups;
    slong * term;
    slong n, k, j, i, l, t;

    tups = flint_malloc(nb * g * sizeof(slong));
    term = flint_malloc(g * sizeof(slong));
    aux = _acb_vec_init(nb);
    acb_init(x);
    fmpz_init(m);
    fmpz_init(p);

    acb_theta_jet_tuples(tups, ord, g);
    for (k = 0; k < nb; k++)
    {
        n = acb_theta_jet_total_order(tups + k * g, g);
        for (j = 0; j < n_pow(g, n); j++)
        {
            for (i = 0; i < g; i++)
            {
                term[i] = 0;
            }
            for (i = 0; i < n; i++)
            {
                term[(j / n_pow(g, i)) % g]++;
            }

            /* multiply by factorials */
            fmpz_one(p);
            for (i = 0; i < g; i++)
            {
                fmpz_fac_ui(m, term[i]);
                fmpz_mul(p, p, m);
            }
            acb_mul_fmpz(x, &v[acb_theta_jet_index(term, g)], p, prec);

            /* view tup as a collection of n indices, enumerate them */
            i = 0;
            for (l = 0; l < g; l++)
            {
                for (t = 0; t < tups[k * g + l]; t++)
                {
                    acb_mul(x, x, acb_mat_entry(N, (j / n_pow(g, i) % g), l), prec);
                    i++;
                }
            }
            acb_add(&aux[k], &aux[k], x, prec);
        }

        fmpz_one(p);
        for (i = 0; i < g; i++)
        {
            fmpz_fac_ui(m, tups[k * g + i]);
            fmpz_mul(p, p, m);
        }
        acb_div_fmpz(&aux[k], &aux[k], p, prec);
    }

    _acb_vec_set(res, aux, nb);

    flint_free(tups);
    flint_free(term);
    _acb_vec_clear(aux, nb);
    acb_clear(x);
    fmpz_clear(m);
    fmpz_clear(p);
}

