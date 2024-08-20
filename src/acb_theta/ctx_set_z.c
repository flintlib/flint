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
acb_theta_ctx_set_z(acb_theta_ctx_t ctx, acb_srcptr z, slong j, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    arb_t u;
    arb_ptr y, t, a;
    acb_ptr s, new_z;
    slong k;

    FLINT_ASSERT(j >= 0 && j < acb_theta_ctx_nb(ctx));

    arb_init(u);
    y = _arb_vec_init(g);
    t = _arb_vec_init(g);
    a = _arb_vec_init(g);
    new_z = _acb_vec_init(g);
    s = _acb_vec_init(g);

    /* We want to compute:
       a - integer vector with even entries which was used for the reduction
       new_z - new values of z after reduction
       c - multiplicative factor for theta values
       v - center of ellipsoid after reduction (if g >= 2)
       as well as the exponentials of the entries of new_z. */

    /* Round t = Yinv y to nearest vector a = 0 mod 2 */
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(t, acb_theta_ctx_yinv(ctx), y, prec);

    _arb_vec_scalar_mul_2exp_si(t, t, g, -1);
    acb_theta_ctx_round(a, t, g);
    _arb_vec_scalar_mul_2exp_si(a, a, g, 1);
    _arb_vec_scalar_mul_2exp_si(t, t, g, 1);
    _arb_vec_set(acb_theta_ctx_as(ctx) + j * g, a, g);

    /* u is exp(pi y^T Yinv y) */
    arb_dot(u, NULL, 0, y, 1, t, 1, g, prec);
    arb_const_pi(acb_theta_ctx_us(ctx) + j, prec);
    arb_mul(acb_theta_ctx_us(ctx) + j, acb_theta_ctx_us(ctx) + j, u, prec);
    arb_exp(acb_theta_ctx_us(ctx) + j, acb_theta_ctx_us(ctx) + j, prec);

    /* v is C(t - a) */
    if (g > 1)
    {
        _arb_vec_sub(t, t, a, g, prec);
        arb_mat_vector_mul_col(acb_theta_ctx_vs(ctx) + j * g, acb_theta_ctx_cho(ctx), t, prec);
    }

    /* new_z is z - tau * a; we further reduce its x-coordinate mod 4 */
    _arb_vec_zero(t, g);
    _acb_vec_set_real_imag(new_z, a, t, g);
    acb_mat_vector_mul_col(new_z, acb_theta_ctx_tau(ctx), new_z, prec);
    _acb_vec_sub(new_z, z, new_z, g, prec);

    _acb_vec_get_real(t, new_z, g);
    _arb_vec_scalar_mul_2exp_si(t, t, g, -2);
    acb_theta_ctx_round(t, t, g);
    _arb_vec_scalar_mul_2exp_si(t, t, g, 2);
    for (k = 0; k < g; k++)
    {
        acb_sub_arb(&new_z[k], &new_z[k], &t[k], prec);
    }

    /* Set exp_z, exp_z_inv, exp_2z, exp_2z_inv */
    for (k = 0; k < g; k++)
    {
        acb_exp_pi_i(&acb_theta_ctx_exp_zs(ctx)[j * g + k], &new_z[k], prec);
        acb_sqr(&acb_theta_ctx_exp_2zs(ctx)[j * g + k], &acb_theta_ctx_exp_zs(ctx)[j * g + k], prec);
        if (acb_is_real(&new_z[k]))
        {
            acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[j * g + k],
                &acb_theta_ctx_exp_zs(ctx)[j * g + k]);
            acb_conj(&acb_theta_ctx_exp_2zs_inv(ctx)[j * g + k],
                &acb_theta_ctx_exp_2zs(ctx)[j * g + k]);
        }
        else
        {
            acb_inv(&acb_theta_ctx_exp_zs_inv(ctx)[j * g + k],
                &acb_theta_ctx_exp_zs(ctx)[j * g + k], prec);
            acb_sqr(&acb_theta_ctx_exp_2zs_inv(ctx)[j * g + k],
                &acb_theta_ctx_exp_zs_inv(ctx)[j * g + k], prec);
        }
    }

    /* c is exp(- i pi a^T (z + new_z)); use new_z as temp */
    _acb_vec_add(new_z, new_z, z, g, prec);
    _arb_vec_zero(t, g);
    _acb_vec_set_real_imag(s, a, t, g);
    acb_dot(&acb_theta_ctx_cs(ctx)[j], NULL, 1, s, 1, new_z, 1, g, prec);
    acb_exp_pi_i(&acb_theta_ctx_cs(ctx)[j], &acb_theta_ctx_cs(ctx)[j], prec);

    arb_clear(u);
    _arb_vec_clear(y, g);
    _arb_vec_clear(t, g);
    _arb_vec_clear(a, g);
    _acb_vec_clear(s, g);
    _acb_vec_clear(new_z, g);
}
