/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_naive_term(acb_t res, acb_srcptr z, const acb_mat_t tau,
    const slong * tup, const slong * n, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr v, w;
    acb_t dot;
    fmpz_t m, t;
    slong k;

    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    acb_init(dot);
    fmpz_init(m);
    fmpz_init(t);

    for (k = 0; k < g; k++)
    {
        acb_set_si(&v[k], n[k]);
    }

    acb_mat_vector_mul_col(w, tau, v, prec);
    acb_dot(res, NULL, 0, v, 1, w, 1, g, prec);
    acb_dot(dot, NULL, 0, v, 1, z, 1, g, prec);
    acb_mul_2exp_si(dot, dot, 1);
    acb_add(res, res, dot, prec);
    acb_exp_pi_i(res, res, prec);

    if (tup != NULL)
    {
        fmpz_one(m);
        for (k = 0; k < g; k++)
        {
            fmpz_set_si(t, n[k]);
            fmpz_pow_ui(t, t, tup[k]);
            fmpz_mul(m, m, t);
        }
        acb_mul_fmpz(res, res, m, prec);
    }

    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    acb_clear(dot);
    fmpz_clear(m);
    fmpz_clear(t);
}
