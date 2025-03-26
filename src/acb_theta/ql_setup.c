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

/* Helper functions */

static void
acb_theta_ctx_tau_copy(acb_theta_ctx_tau_t res, const acb_theta_ctx_tau_t ctx)
{
    slong g = ctx->g;
    slong n = 1 << g;

    FLINT_ASSERT(res->g == g);

    arb_mat_set(&res->yinv, &ctx->yinv);
    arb_mat_set(&res->cho, &ctx->cho);
    acb_mat_set(res->exp_tau_div_4, ctx->exp_tau_div_4);
    acb_mat_set(res->exp_tau_div_2, ctx->exp_tau_div_2);
    acb_mat_set(res->exp_tau, ctx->exp_tau);
    acb_mat_set(res->exp_tau_div_4_inv, ctx->exp_tau_div_4_inv);
    acb_mat_set(res->exp_tau_div_2_inv, ctx->exp_tau_div_2_inv);
    acb_mat_set(res->exp_tau_inv, ctx->exp_tau_inv);

    if (ctx->allow_shift) /* will always be the case here */
    {
        _acb_vec_set(res->exp_tau_a, ctx->exp_tau_a, n * g);
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

static int
acb_theta_contains_even_zero(acb_srcptr vec, slong g)
{
    slong n = 1 << (2 * g);
    slong j;
    int res = 0;

    for (j = 0; j < n; j++)
    {
        if (acb_theta_char_is_even(j, g) && acb_contains_zero(&vec[j]))
        {
            res = 1;
            break;
        }
    }

    return res;
}

/* Find out for which z we can use t=0 at this precision. */
/* Input/output:
   - rts, rts_all, all, nb_steps, nb, distances are as in acb_theta_ql_setup
   - prec should be small
   - is_zero[j] is true iff the corresponding z is zero
   - easy_steps[j] should contain the number of steps that are already known to
     be easy for a given z. This number may increase in the output
   - &vec[j] and ctx_tau should be a valid context for the pairs (z, tau). They
     will get duplicated in acb_theta_ql_setup_easy.
   - At exit, &vec[j] contains the context for 2^k z, where k = easy_steps[j],
     except if easy_steps[j] is nb_steps in which case we guarantee nothing. */

static void
acb_theta_ql_setup_easy(acb_ptr rts, acb_ptr rts_all, slong * easy_steps,
    acb_theta_ctx_z_struct * vec, slong nb, const int * is_zero, acb_theta_ctx_tau_t ctx_tau,
    arb_srcptr distances, slong nb_steps, int all, slong prec)
{
    slong g = ctx_tau->g;
    slong n = 1 << g;
    arb_ptr d;
    acb_ptr th;
    slong j, k;
    int easy;

    d = _arb_vec_init(n);

    for (k = 0; k < nb_steps; k++)
    {
        for (j = 0; j < nb; j++)
        {
            if (easy_steps[j] == nb_steps)
            {
                /* nothing more to be done for this z: we continue without
                   duplicating the context */
                continue;
            }
            else if (easy_steps[j] < k)
            {
                /* steps for this z aren't easy anymore: we continue without
                   duplicating the context */
                continue;
            }
            else if (easy_steps[j] > k)
            {
                /* easy step was set up during a previous pass: we continue but
                   need to update the context */
                acb_theta_ctx_z_dupl(&vec[j], prec + ACB_THETA_LOW_PREC);
                continue;
            }
            /* we now have k == easy_steps[j], and vec[j] contains a valid
               context for 2^k z_j */

            _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, k);
            if (k == 0 && all)
            {
                th = rts_all + j * n * n;
                acb_theta_sum(th, &vec[j], 1, ctx_tau, d, 1, 1, 1, prec);
                /* odd theta constants are known to be zero */
                easy = (is_zero[j] && !acb_theta_contains_even_zero(th, g))
                    || !_acb_vec_contains_zero(th, n * n);
            }
            else
            {
                th = rts + j * (3 * n * nb_steps) + k * (3 * n);
                acb_theta_sum(th, &vec[j], 1, ctx_tau, d, 1, 0, 1, prec);
                easy = !_acb_vec_contains_zero(th, n);
            }

            if (easy)
            {
                /* update easy_steps and duplicate context for k + 1 */
                easy_steps[j]++;
                if (k < nb_steps - 1)
                {
                    acb_theta_ctx_z_dupl(&vec[j], prec + ACB_THETA_LOW_PREC);
                }
            }
        }
        if (k < nb_steps - 1)
        {
            /* duplication on ctx_tau */
            acb_theta_ctx_tau_dupl(ctx_tau, prec + ACB_THETA_LOW_PREC);
        }
    }

    _arb_vec_clear(d, n);
}

/* Find out if the given vector t is usable for hard steps */
/* Input/output:
   - rts, rts_all, all, nb_steps, nb, distances should be as in acb_theta_ql_setup
   - prec should be small
   - t should be an exact real vector of length g
   - easy_steps is as output by acb_theta_ql_setup_easy
   - ctx_tau should contain a valid context for tau; it will get duplicated
   - initially, &vec[j] should contain a valid context for 2^k z where
     k = easy_steps[j]. It will get duplicated, and we guarantee nothing on the
     output.
*/

static int
acb_theta_ql_setup_hard(acb_ptr rts, acb_ptr rts_all, acb_ptr t,
    acb_theta_ctx_z_struct * vec, slong nb, acb_theta_ctx_tau_t ctx_tau,
    const slong * easy_steps, arb_srcptr distances, slong nb_steps, int all, slong prec)
{
    slong g = ctx_tau->g;
    slong n = 1 << g;
    flint_rand_t state;
    acb_theta_ctx_z_struct * aux;
    acb_theta_ctx_z_t ctxt;
    acb_ptr th;
    arb_ptr d;
    slong j, k;
    int res = 1;

    flint_rand_init(state);
    aux = acb_theta_ctx_z_vec_init(2, g);
    acb_theta_ctx_z_init(ctxt, g);
    d = _arb_vec_init(n);

    /* Choose a random t and set context */
    for (k = 0; k < g; k++)
    {
        arb_urandom(acb_realref(&t[k]), state, prec);
        acb_get_mid(&t[k], &t[k]);
    }
    acb_theta_ctx_z_set(ctxt, t, ctx_tau, prec + ACB_THETA_LOW_PREC);

    for (k = 0; (k < nb_steps) && res; k++)
    {
        for (j = 0; (j < nb) && res; j++)
        {
            if (easy_steps[j] > k)
            {
                /* nothing to be done yet for this z: we continue without
                   duplicating the context */
                continue;
            }
            /* we now have k >= easy_steps[j], and vec[j] currently
               contains a valid context for 2^k z_j, while ctxt
               contains a valid context for 2^k t */

            /* set context vector aux at z + t and z + 2t */
            _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, k);
            acb_theta_ctx_z_add_real(&aux[0], &vec[j], ctxt, prec + ACB_THETA_LOW_PREC);
            acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, prec + ACB_THETA_LOW_PREC);

            if (k == 0 && all)
            {
                /* We just need roots for z + 2t */
                th = rts_all + j * n * n;
                acb_theta_sum(th, &aux[1], 1, ctx_tau, d, 1, 1, 1, prec);
                res = !_acb_vec_contains_zero(th, n * n);
            }
            else if (k == 0)
            {
                /* We just need roots for z + 2t */
                th = rts + j * (3 * n * nb_steps) + 2 * n;
                acb_theta_sum(th, &aux[1], 1, ctx_tau, d, 1, 0, 1, prec);
                res = !_acb_vec_contains_zero(th, n);
            }
            else
            {
                /* we need roots at z + t and z + 2t */
                th = rts + j * (3 * n * nb_steps) + k * (3 * n) + n;
                acb_theta_sum(th, aux, 2, ctx_tau, d, 1, 0, 1, prec);
                res = !_acb_vec_contains_zero(th, 2 * n);
            }

            if (k < nb_steps - 1)
            {
                acb_theta_ctx_z_dupl(&vec[j], prec + ACB_THETA_LOW_PREC);
            }
        }
        if (k < nb_steps - 1)
        {
            acb_theta_ctx_tau_dupl(ctx_tau, prec + ACB_THETA_LOW_PREC);
            acb_theta_ctx_z_dupl(ctxt, prec + ACB_THETA_LOW_PREC);
        }
    }

    flint_rand_clear(state);
    acb_theta_ctx_z_vec_clear(aux, 2);
    acb_theta_ctx_z_clear(ctxt);
    _arb_vec_clear(d, n);
    return res;
}

static slong
acb_theta_ql_setup_lost_bits(acb_srcptr rts, acb_srcptr rts_all,
    const slong * easy_steps, const int * is_zero, slong nb,
    arb_srcptr distances, slong nb_steps, int all, slong g)
{
    slong n = 1 << g;
    slong lost_bits, total_lost_bits;
    arb_t x;
    arf_t y, z;
    fmpz_t e;
    slong lp = ACB_THETA_LOW_PREC + nb_steps;
    slong j, k, a, b;

    arb_init(x);
    arf_init(y);
    arf_init(z);
    fmpz_init(e);

    total_lost_bits = 0;
    for (k = 0; k < nb_steps; k++)
    {
        lost_bits = 0;
        for (j = 0; j < nb; j++)
        {
            for (a = 0; a < n; a++)
            {
                arb_zero(x);
                arb_get_lbound_arf(arb_midref(x), &distances[j * n + a], lp);
                arb_mul_2exp_si(x, x, k);
                arb_exp(x, x, lp);
                /* Set y to the minimum absolute value of the roots */
                if (k == 0 && all && easy_steps[j] > 0 && is_zero[j])
                {
                    arf_pos_inf(y);
                    for (b = 0; b < n; b++)
                    {
                        if (acb_theta_char_is_even(n * a + b, g))
                        {
                            acb_get_abs_lbound_arf(z, &rts_all[j * n * n + a * n + b], lp);
                            arf_min(y, y, z);
                        }
                    }
                }
                else if (k == 0 && all) /* easy or hard step, all theta values */
                {
                    arf_pos_inf(y);
                    for (b = 0; b < n; b++)
                    {
                        acb_get_abs_lbound_arf(z, &rts_all[j * n * n + a * n + b], lp);
                        arf_min(y, y, z);
                    }
                }
                else if (k < easy_steps[j]) /* a generic easy step */
                {
                    acb_get_abs_lbound_arf(y, &rts[j * (3 * n * nb_steps) + k * (3 * n) + a], lp);
                }
                else /* a generic hard step */
                {
                    acb_get_abs_lbound_arf(y, &rts[j * (3 * n * nb_steps) + k * (3 * n) + n + a], lp);
                    acb_get_abs_lbound_arf(z, &rts[j * (3 * n * nb_steps) + k * (3 * n) + 2 * n + a], lp);
                    arf_min(y, y, z);
                }
                /* Find out how many bits of precision we approximately lose */
                arb_mul_arf(x, x, y, lp);
                arb_get_lbound_arf(y, x, lp);
                arf_frexp(y, e, y);
                if (fmpz_fits_si(e))
                {
                    lost_bits = FLINT_MAX(lost_bits, - fmpz_get_si(e));
                }
            }
        }
        total_lost_bits += lost_bits + g + 4;
    }

    arb_clear(x);
    arf_clear(y);
    arf_clear(z);
    fmpz_clear(e);
    return total_lost_bits;
}

/* Assume that zs always starts with zero. */

int
acb_theta_ql_setup(acb_ptr rts, acb_ptr rts_all, acb_ptr t, slong * guard, slong * easy_steps,
    acb_srcptr zs, slong nb, const acb_mat_t tau, arb_srcptr distances,
    slong nb_steps, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_ctx_tau_t ctx_tau_1, ctx_tau_2;
    acb_theta_ctx_z_struct * vec;
    slong lowprec, j;
    int done = 0;
    int * is_zero;

    FLINT_ASSERT(nb >= 1);
    FLINT_ASSERT(_acb_vec_is_zero(zs, g));

    for (j = 0; j < nb; j++)
    {
        easy_steps[j] = 0;
    }
    _acb_vec_zero(t, g);
    if (nb_steps == 0)
    {
        *guard = 0;
        return 1;
    }

    acb_theta_ctx_tau_init(ctx_tau_1, 1, g);
    acb_theta_ctx_tau_init(ctx_tau_2, 1, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);
    is_zero = flint_malloc(nb * sizeof(int));

    for (j = 0; j < nb; j++)
    {
        is_zero[j] = _acb_vec_is_zero(zs + j * g, g);
    }

    for (lowprec = 8; (lowprec < prec) && !done; lowprec *= 2)
    {
        /* Set contexts at low precision, but with some additional guard bits */
        acb_theta_ctx_tau_set(ctx_tau_1, tau, lowprec + ACB_THETA_LOW_PREC);
        acb_theta_ctx_tau_copy(ctx_tau_2, ctx_tau_1);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau_1, lowprec + ACB_THETA_LOW_PREC);
        }

        /* Add possible easy steps compared to what we already know */
        acb_theta_ql_setup_easy(rts, rts_all, easy_steps, vec, nb,
            is_zero, ctx_tau_1, distances, nb_steps, all, lowprec);

        /* Let easy_steps[0] be the minimal value (necessary for duplication
           steps) */
        for (j = 1; j < nb; j++)
        {
            easy_steps[0] = FLINT_MIN(easy_steps[0], easy_steps[j]);
        }
        done = (easy_steps[0] == nb_steps);

        /* At this point, for every j such that easy_steps[j] = nb_steps, we
           have computed all the roots we want. Pick an auxiliary t for the
           other hard steps; if it doesn't work, then we restart with an
           increased lowprec. */
        if (!done)
        {
            /* Reset ctx_tau_dupl; note that the input vec in
               setup_hard is exactly as output by setup_easy */
            done = acb_theta_ql_setup_hard(rts, rts_all, t, vec, nb,
                ctx_tau_2, easy_steps, distances, nb_steps, all, lowprec);
        }
    }

    if (done)
    {
        /* Estimate the number of bits lost in the duplication formulas */
        *guard = acb_theta_ql_setup_lost_bits(rts, rts_all, easy_steps,
            is_zero, nb, distances, nb_steps, all, g);
    }

    acb_theta_ctx_tau_clear(ctx_tau_1);
    acb_theta_ctx_tau_clear(ctx_tau_2);
    acb_theta_ctx_z_vec_clear(vec, nb);
    flint_free(is_zero);
    return done;
}
