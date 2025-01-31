/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_QL_TRY 4

static void
acb_theta_ctx_tau_copy(acb_theta_ctx_tau_t res, const acb_theta_ctx_tau_t ctx)
{
    slong g = ctx->g;
    slong n = 1 << g;

    FLINT_ASSERT(res->g == g);

    arb_mat_set(&res->yinv, &ctx->yinv);
    acb_mat_set(res->exp_tau_div_4, ctx->exp_tau_div_4);
    acb_mat_set(res->exp_tau_div_2, ctx->exp_tau_div_2);
    acb_mat_set(res->exp_tau, ctx->exp_tau);

    if (g > 1)
    {
        arb_mat_set(&res->cho, &ctx->cho);
        acb_mat_set(res->exp_tau_div_4_inv, ctx->exp_tau_div_4_inv);
        acb_mat_set(res->exp_tau_div_2_inv, ctx->exp_tau_div_2_inv);
        acb_mat_set(res->exp_tau_inv, ctx->exp_tau_inv);
        _acb_vec_set(res->exp_tau_a_div_2, ctx->exp_tau_a_div_2, n * g);
        _acb_vec_set(res->exp_tau_a, ctx->exp_tau_a, n * g);
        _acb_vec_set(res->exp_tau_a_div_2_inv, ctx->exp_tau_a_div_2_inv, n * g);
        _acb_vec_set(res->exp_tau_a_inv, ctx->exp_tau_a_inv, n * g);
        _acb_vec_set(res->exp_a_tau_a_div_4, ctx->exp_a_tau_a_div_4, n);
    }
}

static int
_acb_vec_contains_zero(acb_srcptr vec, slong nb)
{
    slong k;

    for (k = 0; k < nb; k++)
    {
        if (acb_contains_zero(&vec[k]))
        {
            return 1;
        }
    }
    return 0;
}

/* Assume that zs always starts with zero. */

int
acb_theta_ql_setup(acb_ptr rts, acb_ptr rts_all, acb_ptr t, slong * guard, slong * easy_steps,
    acb_srcptr zs, slong nb, const acb_mat_t tau, arb_srcptr distances,
    slong nb_steps, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau, ctx_tau_dupl;
    acb_theta_ctx_z_t ctxt;
    acb_theta_ctx_z_struct * vec;
    acb_theta_ctx_z_struct * aux;
    flint_rand_t state;
    arb_ptr d;
    slong lowprec;
    slong j, k, l;
    int res, easy;
    int * is_zero;

    FLINT_ASSERT(nb >= 1);
    FLINT_ASSERT(_acb_vec_is_zero(zs, g));

    if (nb_steps == 0)
    {
        *guard = 0;
        for (j = 0; j < nb; j++)
        {
            easy_steps[j] = 0;
        }
        return 1;
    }

    acb_theta_ctx_tau_init(ctx_tau, 1, g);
    acb_theta_ctx_tau_init(ctx_tau_dupl, 1, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);
    aux = acb_theta_ctx_z_vec_init(2, g);
    acb_theta_ctx_z_init(ctxt, g);
    d = _arb_vec_init(n);
    is_zero = flint_malloc(nb * sizeof(int));
    flint_rand_init(state);

    for (j = 0; j < nb; j++)
    {
        is_zero[j] = _acb_vec_is_zero(zs + j * g, g);
    }

    res = 0;
    for (lowprec = 16 + 2 * g; (lowprec < prec) && !res; lowprec *= 2)
    {
        *guard = lowprec;
        /* Set context at precision lowprec */
        /* Find out for which z we can use t=0 at this precision. */
        acb_theta_ctx_tau_set(ctx_tau, tau, lowprec);
        acb_theta_ctx_tau_copy(ctx_tau_dupl, ctx_tau);
        res = 1;
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, lowprec);
            easy_steps[j] = 0;
        }

        for (k = 0; k < nb_steps; k++)
        {
            for (j = 0; j < nb; j++)
            {
                if (easy_steps[j] > k)
                {
                    continue; /* was already set up during a previous pass */
                }
                if (easy_steps[j] < k)
                {
                    continue; /* z was already discarded */
                }

                _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, k);
                if (k == 0 && all)
                {
                    acb_theta_sum(rts_all + j * n * n, &vec[j], 1,
                        ctx_tau_dupl, d, 1, 1, 1, lowprec);
                    /* odd theta constants are known to be zero */
                    if (is_zero[j])
                    {
                        easy = 1;
                        for (l = 0; l < n * n; l++)
                        {
                            if (acb_theta_char_is_even(l, g))
                            {
                                if(acb_contains_zero(&rts_all[j * n * n + l]))
                                {
                                    easy = 0;
                                    break;
                                }
                            }
                            else
                            {
                                acb_zero(&rts_all[j * n * n + l]);
                            }
                        }
                    }
                    else
                    {
                        easy = !_acb_vec_contains_zero(rts_all + j * n * n, n * n);
                    }
                }
                else
                {
                    acb_theta_sum(rts + j * (3 * n * nb_steps) + k * (3 * n),
                        &vec[j], 1, ctx_tau_dupl, d, 1, 0, 1, lowprec);
                    easy = !_acb_vec_contains_zero(rts + j * (3 * n * nb_steps) + k * (3 * n), n);
                }
                if (easy)
                {
                    easy_steps[j]++;
                    if (k < nb_steps - 1)
                    {
                        acb_theta_ctx_z_dupl(&vec[j], lowprec);
                    }
                }
                else
                {
                    res = 0; /* trigger search for t */
                }
            }
            if (k < nb_steps - 1)
            {
                acb_theta_ctx_tau_dupl(ctx_tau_dupl, lowprec);
            }
        }

        /* At this point, for every j such that easy_steps[j] = nb_steps,
           we have computed all the roots we want. Pick an auxiliary t for the
           other j's; we must include 0 in the "not done" list */
        _acb_vec_zero(t, g);
        for (j = 0; j < nb; j++)
        {
            easy_steps[0] = FLINT_MIN(easy_steps[0], easy_steps[j]);
        }

        for (l = 0; (l < ACB_THETA_QL_TRY) && !res; l++)
        {
            res = 1;
            for (k = 0; k < g; k++)
            {
                arb_urandom(acb_realref(&t[k]), state, prec);
                acb_get_mid(&t[k], &t[k]);
            }
            /* Reinitialize contexts */
            acb_theta_ctx_tau_copy(ctx_tau_dupl, ctx_tau);
            acb_theta_ctx_z_set(ctxt, t, ctx_tau, lowprec);

            /* Find out if roots are all nonzero */
            for (k = 0; (k < nb_steps) && res; k++)
            {
                for (j = 0; (j < nb) && res; j++)
                {
                    if (k < easy_steps[j])
                    {
                        continue; /* this step was already handled. */
                    }

                    acb_theta_ctx_z_add_real(&aux[0], &vec[j], ctxt, lowprec);
                    acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, lowprec);
                    _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, k);
                    if (k == 0 && all)
                    {
                        /* We just need roots for z + 2t */
                        acb_theta_sum(rts_all + j * n * n, aux + 1,
                            1, ctx_tau_dupl, d, 1, 1, 1, lowprec);
                        res = res && !_acb_vec_contains_zero(rts_all + j * n * n, n * n);

                        /* if (!res)
                        {
                            flint_printf("(ql_setup) fail at guard = %wd (l = %wd, k = %wd, j = %wd)\n", lowprec, l, k, j);
                            _acb_vec_printd(rts_all + j * n * j, n * n, 5);
                            }*/
                    }
                    else
                    {
                        acb_theta_sum(rts + j * (3 * n * nb_steps) + k * (3 * n) + n,
                            aux, 2, ctx_tau_dupl, d, 1, 0, 1, lowprec);
                        res = res && !_acb_vec_contains_zero(rts + j * (3 * n * nb_steps)
                            + k * (3 * n) + n, 2 * n);

                        /*if (!res)
                        {
                            flint_printf("(ql_setup) fail at guard = %wd (l = %wd, k = %wd, j = %wd)\n", lowprec, l, k, j);
                            _acb_vec_printd(rts + j * (3 * n * nb_steps) + k * (3 * n) + n, 2 * n, 5);
                            }*/
                    }
                    if (k < nb_steps - 1)
                    {
                        acb_theta_ctx_z_dupl(&vec[j], lowprec);
                    }
                }
                if (k < nb_steps - 1)
                {
                    acb_theta_ctx_tau_dupl(ctx_tau_dupl, lowprec);
                    acb_theta_ctx_z_dupl(ctxt, lowprec);
                }
            }
            /* If res, then all the roots are computed, and we are done.
               Otherwise, try a new vector t and reinitialize the entries of vec. */
            if (!res)
            {
                for (j = 0; j < nb; j++)
                {
                    if (easy_steps[j] == nb_steps)
                    {
                        continue;
                    }
                    acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, lowprec);
                    for (k = 0; k < easy_steps[j]; k++)
                    {
                        acb_theta_ctx_z_dupl(&vec[j], lowprec);
                    }
                }
            }
        }
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_tau_clear(ctx_tau_dupl);
    acb_theta_ctx_z_vec_clear(vec, nb);
    acb_theta_ctx_z_vec_clear(aux, 2);
    acb_theta_ctx_z_clear(ctxt);
    _arb_vec_clear(d, n);
    flint_free(is_zero);
    flint_rand_clear(state);
    return res;
}
