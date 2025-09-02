/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb.h"
#include "arb_mat.h"

/* Adapted from fmpz_mat_is_reduced_gram */
int
arb_mat_spd_is_lll_reduced(const arb_mat_t A, double delta, double eta, slong prec)
{
    slong d = arb_mat_nrows(A);
    arb_mat_t r, mu;
    arb_ptr s;
    arb_t delta_arb, eta_arb, t;
    slong i, j, k;
    int res = 1;

    if (d <= 1)
    {
        return 1;
    }

    arb_mat_init(r, d, d);
    arb_mat_init(mu, d, d);
    s = _arb_vec_init(d);
    arb_init(delta_arb);
    arb_init(eta_arb);
    arb_init(t);

    arb_set_d(delta_arb, delta);
    arb_set_d(eta_arb, eta);

    arb_set(arb_mat_entry(r, 0, 0), arb_mat_entry(A, 0, 0));

    for (i = 1; (i < d) && res; i++)
    {
        arb_set(&s[0], arb_mat_entry(A, i, i));
        for (j = 0; (j < i) && res; j++)
        {
            arb_set(arb_mat_entry(r, i, j), arb_mat_entry(A, i, j));
            for (k = 0; k < j; k++)
            {
                arb_submul(arb_mat_entry(r, i, j), arb_mat_entry(mu, j, k),
                    arb_mat_entry(r, i, k), prec);
            }
            arb_div(arb_mat_entry(mu, i, j), arb_mat_entry(r, i, j),
                arb_mat_entry(r, j, j), prec);
            arb_abs(t, arb_mat_entry(mu, i, j));
            if (!arb_le(t, eta_arb))
            {
                res = 0;
            }
            arb_set(&s[j + 1], &s[j]);
            arb_submul(&s[j + 1], arb_mat_entry(mu, i, j), arb_mat_entry(r, i, j), prec);
        }
        arb_set(arb_mat_entry(r, i, i), &s[i]);
        arb_mul(t, delta_arb, arb_mat_entry(r, i - 1, i - 1), prec);
        if (!arb_le(t, &s[i - 1]))
        {
            res = 0;
        }
    }

    arb_mat_clear(r);
    arb_mat_clear(mu);
    _arb_vec_clear(s, d);
    arb_clear(delta_arb);
    arb_clear(eta_arb);
    arb_clear(t);
    return res;
}
