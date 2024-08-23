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

static void
acb_theta_ql_step_1(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;

    acb_theta_agm_mul_tight(res, th0, th, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(res, res, n, g);
    acb_theta_agm_sqrt(res, res, rts, n, prec);
}

static void
acb_theta_ql_step_2(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;

    aux = _acb_vec_init(3 * n);

    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    _acb_vec_set(res, aux, 3 * n);
    _acb_vec_clear(aux, 3 * n);
}


static void
acb_theta_ql_step_3(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;
    ulong a;

    aux = _acb_vec_init(3 * n);

    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    acb_theta_agm_mul_tight(aux, th0 + n, th + n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux, aux, n, g);
    for (a = 0; a < n; a++)
    {
        acb_div(&aux[a], &aux[a], &aux[2 * n + a], prec);
    }

    _acb_vec_set(res, aux, 3 * n);
    _acb_vec_clear(aux, 3 * n);
}

static void
acb_theta_ql_one_step(acb_ptr th, acb_srcptr rts, arb_srcptr d0, arb_srcptr d,
    int t_is_zero, int z_is_zero, int last_step, slong g, slong prec)
{
    slong n = 1 << g;
    slong nb = (t_is_zero ? 1 : 3) * (z_is_zero ? 2 : 1);
    acb_ptr next;

    next = _acb_vec_init(nb * n);

    if (t_is_zero)
    {
        acb_theta_ql_step_1(next, th, th, rts, d0, d0, g, prec);
        if (!z_is_zero)
        {
            acb_theta_ql_step_1(next + n, th, th + n,
                rts + n, d0, d, g, prec);
        }
    }
    else /* t not zero */
    {
        acb_theta_ql_step_3(next, th, th, rts, d0, d0, g, prec);
        if (!z_is_zero && last_step)
        {
            acb_theta_ql_step_3(next + 3 * n, th, th + 3 * n,
                rts + 2 * n, d0, d, g, prec);
        }
        else if (!z_is_zero)
        {
            acb_theta_ql_step_2(next + 3 * n, th, th + 3 * n,
                rts + 2 * n, d0, d, g, prec);
        }
    }

    _acb_vec_set(th, next, nb * n);
    _acb_vec_clear(next, nb * n);
}
/*
int
acb_theta_ql_steps(acb_ptr th, slong nb, slong split, const acb_theta_ctx_t ctx,
    slong guard, slong prec)
    {*/
    /* slong g = acb_theta_ctx_g(tau); */
    /* slong n = 1 << g; */
    /* acb_theta_ctx_t new_ctx; */
    /* acb_ptr rts, th; */
    /* arb_ptr d0, d; */
    /* slong nb_rts = (ctx->t_is_zero ? 1 : 2) * (ctx->z_is_zero ? 1 : 2); */
    /* slong nb_th = (ctx->t_is_zero ? 1 : 3) * (ctx->z_is_zero ? 1 : 2); */
    /* slong j, k; */
    /* int res = 1; */

    /* acb_theta_ctx_init(new_ctx, 6, g); */
    /* rts = _acb_vec_init(nb * n * nb_rts); */
    /* th = _acb_vec_init(nb_th); */
    /* d0 = _arb_vec_init(n); */
    /* d = _arb_vec_init(n); */

    /* acb_theta_ctx_copy(new_ctx, ctx); */
    /* for (j = 0; (j < nb) && res; j++) */
    /* { */
    /*     /\* Compute roots for step number j. */
    /*        - if t = 0, we need theta_a0(0, tau) and perhaps theta_a0(z, tau) */
    /*        - if t is not zero, we need theta_a0(t, tau), theta_a0(2t, tau) */
    /*        and perhaps theta_a0(z + t, tau), theta_a0(z + 2t, tau) */
    /*        TODO: when z is real, compute only one ellipsoid *\/ */
    /*     if (ctx->t_is_zero) */
    /*     { */
    /*         if (ctx->z_is_zero) */
    /*         { */
    /*             acb_theta_sum_a0(rts + j * nb_rts * n, new_ctx, 0, 1, 1, guard); */
    /*         } */
    /*         else */
    /*         { */
    /*             /\* todo: only one call for real z? *\/ */
    /*             acb_theta_sum_a0(rts + j * nb_rts * n, new_ctx, 0, 1, 1, guard); */
    /*             acb_theta_sum_a0(rts + j * nb_rts * n + n, new_ctx, 3, 1, ctx->z_is_real, guard); */
    /*         } */
    /*     } */
    /*     else /\* t not zero *\/ */
    /*     { */
    /*         if (ctx->z_is_zero) */
    /*         { */
    /*             acb_theta_sum_a0(rts + j * nb_rts * n, new_ctx, 1, 2, 1, guard); */
    /*         } */
    /*         else */
    /*         { */
    /*             acb_theta_sum_a0(rts + j * nb_rts * n, new_ctx, 1, 2, 1, guard); */
    /*             acb_theta_sum_a0(rts + j * nb_rts * n + 2 * n, new_ctx, 4, 2, ctx->z_is_real, guard); */
    /*         } */
    /*     } */
    /*     /\* Decide whether to continue: roots must not contain zero *\/ */
    /*     for (k = 0; (k < n * nb_rts) && res; k++) */
    /*     { */
    /*         res = res && !acb_contains_zero(&rts[j * nb_rts * n + k]); */
    /*     } */
    /*     /\* Duplicate context at high precision *\/ */
    /*     if (res) */
    /*     { */
    /*         acb_theta_ctx_dupl(new_ctx, prec); */
    /*     } */
    /* } */

    /* if (res && split == 0) */
    /* { */
    /*     /\* Initialize using sum_a0 at high precision *\/ */
    /*     if (ctx->t_is_zero) */
    /*     { */
    /*         /\* TODO: only one call if z is real *\/ */
    /*         if (ctx->z_is_zero) */
    /*         { */
    /*             acb_theta_sum_a0(th, new_ctx, 0, 1, 1, prec); */
    /*         } */
    /*         else */
    /*         { */
    /*             acb_theta_sum_a0(th, new_ctx, 0, 1, 1, prec); */
    /*             acb_theta_sum_a0(th + n, new_ctx, 3, 1, ctx->z_is_real, prec); */
    /*         } */
    /*     } */
    /*     else /\* t not zero *\/ */
    /*     { */
    /*         if (ctx->z_is_real) */
    /*         { */
    /*             acb_theta_sum_a0(th, new_ctx, 0, nb_th, 1, prec); */
    /*         } */
    /*         else */
    /*         { */
    /*             acb_theta_sum_a0(th, new_ctx, 0, 3, 1, prec); */
    /*             acb_theta_sum_a0(th, new_ctx, 3, 3, 0, prec); */
    /*         } */
    /*     } */
    /* } */
    /* else if (res) */
    /* { */
    /*     /\* Initialize with splitting startegy *\/ */
    /*     res = acb_theta_ql_split(th, split, new_ctx, guard, prec); */
    /* } */

    /* /\* Perform duplication steps if successful so far *\/ */
    /* if (res) */
    /* { */
    /*     for (j = nb - 1; j >= 0; j--) */
    /*     { */
    /*         _arb_vec_scalar_mul_2exp_si(d0, acb_theta_ctx_d0(ctx), n, j); */
    /*         _arb_vec_scalar_mul_2exp_si(d, acb_theta_ctx_d(ctx), n, j); */
    /*         acb_theta_ql_one_step(th, rts + j * nb_rts * n, d0, d, */
    /*             ctx->t_is_zero, ctx->z_is_zero, j == 0, g, prec); */
    /*     } */
    /* } */

    /* /\* Do we need this f factor ? */
    /* acb_siegel_yinv(Yinv, tau, prec); */
    /* _acb_vec_get_imag(y, z, g); */
    /* arb_mat_vector_mul_col(w, Yinv, y, prec); */
    /* arb_dot(f, NULL, 1, y, 1, w, 1, g, prec); */
    /* arb_const_pi(c, prec); */
    /* arb_mul(f, f, c, prec); *\/ */

    /* acb_theta_ctx_clear(new_ctx); */
    /* _acb_vec_init(rts, nb * n * nb_rts); */
    /* _acb_vec_init(th, nb_th); */
    /* _arb_vec_init(d0, n); */
    /* _arb_vec_init(d, n); */
    /* return res; */
/*}*/

