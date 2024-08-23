/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_CTX_MAX_VMAG 1000000

static void
acb_theta_ctx_round(arb_ptr a, arb_srcptr v, slong g)
{
    slong j;
    fmpz_t m;

    fmpz_init(m);

    for (j = 0; j < g; j++)
    {
        if (arb_is_finite(&v[j])
            && arf_cmpabs_2exp_si(arb_midref(&v[j]), ACB_THETA_CTX_MAX_VMAG) <= 0)
        {
            arf_get_fmpz(m, arb_midref(&v[j]), ARF_RND_NEAR);
            arb_set_fmpz(&a[j], m);
        }
        else
        {
            arb_zero(&a[j]);
        }
    }

    fmpz_clear(m);
}

void
acb_theta_ctx_z_set(acb_theta_ctx_z_t ctx, acb_srcptr z, const acb_theta_ctx_tau_t ctx_tau, slong prec)
{
    slong g = acb_theta_ctx_g(ctx_tau);
    arb_t u;
    arb_ptr y, t, r;
    acb_ptr s, new_z;
    slong k;
    int is_real;

    arb_init(u);
    y = _arb_vec_init(g);
    t = _arb_vec_init(g);
    r = _arb_vec_init(g);
    new_z = _acb_vec_init(g);
    s = _acb_vec_init(g);

    /* We want to compute:
       r - integer vector with even entries which was used for the reduction
       new_z - new values of z after reduction
       c - multiplicative factor for theta values
       v - center of ellipsoid after reduction (if g >= 2)
       as well as the exponentials of the entries of new_z. */

    /* Round t = Yinv y to nearest vector r = 0 mod 2 */
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(t, acb_theta_ctx_yinv(ctx_tau), y, prec);

    _arb_vec_scalar_mul_2exp_si(t, t, g, -1);
    acb_theta_ctx_round(r, t, g);
    _arb_vec_scalar_mul_2exp_si(r, r, g, 1);
    _arb_vec_scalar_mul_2exp_si(t, t, g, 1);
    _arb_vec_set(acb_theta_ctx_r(ctx), r, g);

    /* u is exp(pi y^T Yinv y) and v is C(t - r) */
    if (g > 1)
    {
        arb_dot(u, NULL, 0, y, 1, t, 1, g, prec);
        arb_const_pi(acb_theta_ctx_u(ctx), prec);
        arb_mul(acb_theta_ctx_u(ctx), acb_theta_ctx_u(ctx), u, prec);
        arb_exp(acb_theta_ctx_u(ctx), acb_theta_ctx_u(ctx), prec);

        _arb_vec_sub(t, t, r, g, prec);
        arb_mat_vector_mul_col(acb_theta_ctx_v(ctx), acb_theta_ctx_cho(ctx_tau), t, prec);
    }

    /* new_z is z - tau * r */
    _arb_vec_zero(t, g);
    _acb_vec_set_real_imag(new_z, r, t, g);
    acb_mat_vector_mul_col(new_z, acb_theta_ctx_tau(ctx_tau), new_z, prec);
    _acb_vec_sub(new_z, z, new_z, g, prec);

    /* Set z_is_real, exp_z, exp_z_inv, exp_2z, exp_2z_inv */
    for (k = 0; k < g; k++)
    {
        acb_exp_pi_i(&acb_theta_ctx_exp_z(ctx)[k], &new_z[k], prec);
        if (g > 1)
        {
            is_real = acb_is_real(&new_z[k]);
            acb_sqr(&acb_theta_ctx_exp_2z(ctx)[k], &acb_theta_ctx_exp_z(ctx)[k], prec);
            acb_theta_ctx_exp_inv(&acb_theta_ctx_exp_z_inv(ctx)[k],
                &acb_theta_ctx_exp_z(ctx)[k], &new_z[k], is_real, prec);
            acb_theta_ctx_sqr_inv(&acb_theta_ctx_exp_2z_inv(ctx)[k],
                &acb_theta_ctx_exp_z_inv(ctx)[k], &acb_theta_ctx_exp_2z(ctx)[k],
                is_real, prec);
        }
    }
    _acb_vec_set(acb_theta_ctx_z(ctx), new_z, g);
    acb_theta_ctx_is_real(ctx) = _acb_vec_is_real(new_z, g);

    /* c is exp(- i pi r^T (z + new_z)); use new_z as temp */
    _acb_vec_add(new_z, new_z, z, g, prec);
    _arb_vec_zero(t, g);
    _acb_vec_set_real_imag(s, r, t, g);
    acb_dot(acb_theta_ctx_c(ctx), NULL, 1, s, 1, new_z, 1, g, prec);
    acb_exp_pi_i(acb_theta_ctx_c(ctx), acb_theta_ctx_c(ctx), prec);

    arb_clear(u);
    _arb_vec_clear(y, g);
    _arb_vec_clear(t, g);
    _arb_vec_clear(r, g);
    _acb_vec_clear(s, g);
    _acb_vec_clear(new_z, g);
}
