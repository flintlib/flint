/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

void acb_theta_one_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, ulong ab, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong j;

    FLINT_ASSERT(nb >= 0);
    FLINT_ASSERT(ab >= 0 && ab < (1 << (2 * g)));

    if (g == 1)
    {
        /* call acb_modular_theta_sum directly */
        acb_t w, q14, q;
        acb_ptr res;

        acb_init(w);
        acb_init(q14);
        acb_init(q);
        res = _acb_vec_init(4);

        acb_mul_2exp_si(q14, acb_mat_entry(tau, 0, 0), -2);
        acb_exp_pi_i(q14, q14, prec);
        acb_pow_ui(q, q14, 4, prec);

        for (j = 0; j < nb; j++)
        {
            acb_exp_pi_i(w, &zs[j], prec);
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                w, acb_is_real(&zs[j]), q, 1, prec);
            if (ab == 0)
            {
                acb_set(&th[j], &res[2]);
            }
            else if (ab == 1)
            {
                acb_set(&th[j], &res[3]);
            }
            else if (ab == 2)
            {
                acb_set(&th[j], &res[1]);
            }
            else
            {
                acb_neg(&th[j], &res[0]);
            }
            if (ab >= 2)
            {
                acb_mul(&th[j], &th[j], q14, prec);
            }
        }

        acb_clear(w);
        acb_clear(q14);
        acb_clear(q);
        _acb_vec_clear(res, 4);
    }
    else if (ab == 0)
    {
        /* Call 00_notransform directly */
        acb_theta_00_notransform(th, zs, nb, tau, prec);
    }
    else
    {
        /* theta_ab(z, tau) = exp(pi i a^T tau a/4) exp(2 pi i a^T (z + b/2))
           theta_00(z + tau a/2 + b/2, tau) */
        acb_ptr new_zs, v, w;
        acb_t c, x;
        ulong b = ab % (1 << g);
        ulong a = ab >> g;

        new_zs = _acb_vec_init(nb * g);
        v = _acb_vec_init(g);
        w = _acb_vec_init(g);
        acb_init(c);
        acb_init(x);

        acb_theta_char_get_acb(v, a, g);
        acb_mat_vector_mul_col(v, tau, v, prec); /* tau.a/2 */
        acb_theta_char_get_acb(w, b, g);
        _acb_vec_add(w, v, w, g, prec);
        for (j = 0; j < nb; j++)
        {
            _acb_vec_add(new_zs + j * g, zs + j * g, w, g, prec);
        }

        acb_theta_00_notransform(th, new_zs, nb, tau, prec);

        acb_theta_char_dot_acb(c, a, v, g, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_char_get_acb(w, b, g);
            _acb_vec_add(w, w, zs + j * g, g, prec);
            acb_theta_char_dot_acb(x, a, w, g, prec);
            acb_mul_2exp_si(x, x, 1);
            acb_add(x, x, c, prec);
            acb_exp_pi_i(x, x, prec);
            acb_mul(&th[j], &th[j], x, prec);
        }

        _acb_vec_clear(new_zs, nb * g);
        _acb_vec_clear(v, g);
        _acb_vec_clear(w, g);
        acb_clear(c);
        acb_clear(x);
    }
}
