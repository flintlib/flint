/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int
acb_theta_round(arb_ptr a, arb_srcptr v, slong g)
{
    slong j;
    fmpz_t m;
    int res = 1;

    fmpz_init(m);

    for (j = 0; (j < g) && res; j++)
    {
        res = arb_is_finite(&v[j])
            && arf_cmpabs_2exp_si(arb_midref(&v[j]), 1000000) <= 0;
        if (res)
        {
            arf_get_fmpz(m, arb_midref(&v[j]), ARF_RND_NEAR);
            arb_set_fmpz(&a[j], m);
        }
    }

    fmpz_clear(m);
    return res;
}

int
acb_theta_reduce_z(acb_ptr new_zs, arb_ptr rs, acb_ptr cs, acb_srcptr zs,
    slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t cho, yinv;
    arb_ptr y;
    acb_ptr t, x;
    slong j, k;
    int res = 1;

    FLINT_ASSERT(nb >= 0);

    arb_mat_init(cho, g, g);
    arb_mat_init(yinv, g, g);
    y = _arb_vec_init(g);
    t = _acb_vec_init(g);
    x = _acb_vec_init(g);

    acb_siegel_cho(cho, tau, prec);
    acb_siegel_yinv(yinv, tau, prec);

    for (j = 0; j < nb; j++)
    {
        /* Round Yinv y to nearest vector r = 0 mod 2 */
        _acb_vec_get_imag(y, zs + j * g, g);
        arb_mat_vector_mul_col(y, yinv, y, prec);
        _arb_vec_scalar_mul_2exp_si(y, y, g, -1);
        res = acb_theta_round(rs + j * g, y, g);
        _arb_vec_scalar_mul_2exp_si(rs + j * g, rs + j * g, g, 1);
        if (!res)
        {
            break;
        }

        /* x = new_z is z - tau * r */
        _arb_vec_zero(y, g);
        _acb_vec_set_real_imag(x, rs + j * g, y, g);
        acb_mat_vector_mul_col(x, tau, x, prec);
        _acb_vec_sub(x, zs + j * g, x, g, prec);

        /* c is exp(- i pi r^T (z + x)) */
        _acb_vec_add(t, x, zs + j * g, g, prec);
        _acb_vec_get_real(y, t, g);
        arb_dot(acb_realref(&cs[j]), NULL, 1, rs + j * g, 1, y, 1, g, prec);
        _acb_vec_get_imag(y, t, g);
        arb_dot(acb_imagref(&cs[j]), NULL, 1, rs + j * g, 1, y, 1, g, prec);
        acb_exp_pi_i(&cs[j], &cs[j], prec);

        /* Further reduce real part of x modulo 2 */
        _acb_vec_get_real(y, x, g);
        _arb_vec_scalar_mul_2exp_si(y, y, g, -1);
        res = acb_theta_round(y, y, g);
        if (res)
        {
            _arb_vec_scalar_mul_2exp_si(y, y, g, 1);
            for (k = 0; k < g; k++)
            {
                acb_sub_arb(&x[k], &x[k], &y[k], prec);
            }
        }
        else
        {
            /* Still OK; set real part to [0,2] */
            for (k = 0; k < g; k++)
            {
                arb_unit_interval(acb_realref(&x[k]));
                arb_mul_2exp_si(acb_realref(&x[k]), acb_realref(&x[k]), 1);
            }
            res = 1;
        }

        /* Set new_z */
        _acb_vec_set(new_zs + j * g, x, g);
    }

    arb_mat_clear(cho);
    arb_mat_clear(yinv);
    _arb_vec_clear(y, g);
    _acb_vec_clear(t, g);
    _acb_vec_clear(x, g);
    return res;
}
