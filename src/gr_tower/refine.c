/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/*
    Given a monic factor g of the minimal polynomial m_k, decides
    numerically which of g and m_k / g vanishes at the enclosure of a_k,
    and replaces m_k by that factor.
*/
/* (the enclosure of the generator a_k, for the factor selection) */
typedef struct { gr_tower_struct * T; slong k; } _get_z_step_arg;

static int
_get_z_step(acb_t z, slong prec, void * arg)
{
    _get_z_step_arg * a = (_get_z_step_arg *) arg;
    return gr_tower_step_get_acb(z, a->T, a->k, prec);
}

int
gr_tower_refine_step(gr_tower_t T, slong k, const gr_poly_t g)
{
    gr_tower_gen_struct * step = GR_TOWER_STEP(T, k - 1);
    gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
    const gr_poly_struct * m = gr_poly_quotient_ctx_modulus(step->ctx);
    gr_poly_t h, r;
    acb_t z;
    int status = GR_SUCCESS;
    int which = -1;

    if (g->length < 2 || g->length >= m->length)
        return GR_DOMAIN;

    gr_poly_init(h, below);
    gr_poly_init(r, below);
    acb_init(z);

    /* h = m / g */
    status |= gr_poly_divrem(h, r, m, g, below);

    if (status == GR_SUCCESS && gr_poly_is_zero(r, below) != T_TRUE)
        status = GR_UNABLE;   /* not actually a factor (should not happen) */

    /* the factor among g, h vanishing at the generator */
    if (status == GR_SUCCESS)
    {
        gr_ctx_t pctx;
        gr_vec_t fac;
        _get_z_step_arg a = { T, k };
        gr_ctx_init_gr_poly(pctx, below);
        gr_vec_init(fac, 2, pctx);
        status |= gr_poly_set(gr_vec_entry_ptr(fac, 0, pctx), g, below);
        gr_poly_swap(gr_vec_entry_ptr(fac, 1, pctx), h, below);
        which = _gr_tower_select_poly_factor(fac, k - 1, _get_z_step, &a, T);
        gr_poly_swap(h, gr_vec_entry_ptr(fac, 1, pctx), below);
        gr_vec_clear(fac, pctx);
        gr_ctx_clear(pctx);
    }

    if (status == GR_SUCCESS && which == -1)
        status = GR_UNABLE;

    if (status == GR_SUCCESS)
    {
        const gr_poly_struct * f = (which == 0) ? g : h;

        status |= gr_poly_quotient_ctx_refine(step->ctx, f);

        if (status == GR_SUCCESS)
        {
            gr_poly_quotient_ctx_clear_zero_divisors(step->ctx);

            /* re-derive the enclosure from the new minimal polynomial:
               for a linear factor, the root is exact in F_{k-1} */
            if (f->length == 2)
            {
                status |= gr_tower_get_acb_at(z, gr_poly_coeff_srcptr(f, 0, below), k - 1, GR_TOWER_DEFAULT_PREC, T);
                acb_neg(z, z);
                /* (both enclose the root: the intersection keeps the
                   isolation of the old enclosure) */
                if (status == GR_SUCCESS && acb_overlaps(z, &step->enclosure))
                {
                    arb_intersection(acb_realref(z), acb_realref(z), acb_realref(&step->enclosure), GR_TOWER_DEFAULT_PREC);
                    arb_intersection(acb_imagref(z), acb_imagref(z), acb_imagref(&step->enclosure), GR_TOWER_DEFAULT_PREC);
                }
                acb_set(&step->enclosure, z);
                step->enclosure_prec = GR_TOWER_DEFAULT_PREC;
            }

            T->version++;
            T->moduli_version++;
        }
    }

    gr_poly_clear(h, below);
    gr_poly_clear(r, below);
    acb_clear(z);

    return status;
}

int
gr_tower_refine(gr_tower_t T)
{
    slong k;
    int changed = 0, progress;
    int rounds = 0;

    do
    {
        progress = 0;

        for (k = 1; k <= T->length; k++)
        {
            gr_ctx_struct * ctx = GR_TOWER_STEP_CTX(T, k - 1);
            gr_ctx_struct * below = gr_tower_field_at(T, k - 1);

            while (gr_poly_quotient_ctx_num_zero_divisors(ctx) > 0)
            {
                gr_poly_t g;
                const gr_poly_struct * m;
                int status;
                slong j, lower_pending = 0;

                gr_poly_init(g, below);
                GR_MUST_SUCCEED(gr_poly_set(g, gr_poly_quotient_ctx_zero_divisor(ctx, 0), below));

                /* the modulus may have been refined since g was recorded */
                m = gr_poly_quotient_ctx_modulus(ctx);
                status = gr_poly_gcd(g, g, m, below);

                if (status == GR_SUCCESS)
                {
                    if (g->length >= 2 && g->length < m->length)
                    {
                        status = gr_tower_refine_step(T, k, g);   /* clears the table on success */
                        if (status == GR_SUCCESS)
                            changed = progress = 1;
                    }

                    if (status != GR_SUCCESS || !(g->length >= 2 && g->length < m->length))
                    {
                        /* unusable entry (trivial gcd or numerical failure): drop it */
                        gr_poly_quotient_ctx_clear_zero_divisors(ctx);
                    }
                }
                else
                {
                    /* the gcd computation itself may have exposed a zero
                       divisor at a lower level; if so, let the next round
                       handle it first, then come back to this level */
                    for (j = 1; j < k; j++)
                        lower_pending += gr_poly_quotient_ctx_num_zero_divisors(GR_TOWER_STEP_CTX(T, j - 1));

                    if (lower_pending)
                        progress = 1;
                    else
                        gr_poly_quotient_ctx_clear_zero_divisors(ctx);

                    gr_poly_clear(g, below);
                    break;
                }

                gr_poly_clear(g, below);
            }
        }

        rounds++;
    }
    while (progress && rounds < 1000);

    return changed;
}

/*
    Zero test in the field F_k (dynamic evaluation), without any
    consideration of conjectural transcendence: returns
    GR_TOWER_ZERO if x is zero, GR_TOWER_NONZERO if x is numerically
    separated from zero, GR_TOWER_FIELD_NONZERO if x is a nonzero element
    of the field (a nonzero rational function of the transcendental
    generators over the algebraic part), and GR_TOWER_UNKNOWN otherwise.
*/
int
_gr_tower_field_is_zero_at(gr_srcptr x, slong k, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field_at(T, k);
    gr_ptr t;
    int res;
    int status;

    GR_TMP_INIT(t, top);

    for (;;)
    {
        if (gr_is_zero(x, top) == T_TRUE)
        {
            res = GR_TOWER_ZERO;
            break;
        }

        /* Cheap numerical separation from zero: sound regardless of
           whether the tower is a field, and avoids computing an inverse
           (which can be enormous even when x is small). */
        {
            acb_t z;
            slong prec;
            int separated = 0;

            /* (a nonzero element is separated at a moderate precision;
               the exact test below, an inversion in the nested field, can
               be far more expensive than many more bits, so the precision
               is pushed well beyond the default first) */
            acb_init(z);
            for (prec = GR_TOWER_DEFAULT_PREC; prec <= GR_TOWER_OPTION(T, GR_TOWER_OPT_PREC_LIMIT); prec *= 2)
            {
                if (gr_tower_get_acb_at(z, x, k, prec, T) == GR_SUCCESS && !acb_contains_zero(z))
                {
                    separated = 1;
                    break;
                }
            }
            acb_clear(z);

            if (separated)
            {
                res = GR_TOWER_NONZERO;
                break;
            }
        }

        if (k == 0)
        {
            res = GR_TOWER_FIELD_NONZERO;
            break;
        }

        /* with all moduli up to level k proved irreducible, the chain is
           a field and a nonzero element is invertible: no need to
           compute the inverse (which can be enormous) */
        {
            slong j;
            int all_proven = 1;
            for (j = 1; j <= k && all_proven; j++)
                if (GR_TOWER_STEP(T, j - 1)->status != GR_TOWER_STATUS_PROVEN)
                    all_proven = 0;
            if (!all_proven)
            {
                /* a modular proof is far cheaper than an inversion (and
                   Trager's method may refine the tower instead: then the
                   element is tested again modulo the new moduli) */
                ulong version = T->version;
                gr_tower_prove_modular(T, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES));
                if (T->version != version)
                    continue;
                all_proven = 1;
                for (j = 1; j <= k && all_proven; j++)
                    if (GR_TOWER_STEP(T, j - 1)->status != GR_TOWER_STATUS_PROVEN)
                        all_proven = 0;
            }
            if (all_proven)
            {
                res = GR_TOWER_FIELD_NONZERO;
                break;
            }
        }

        status = gr_inv(t, x, top);

        if (status == GR_SUCCESS)
        {
            res = GR_TOWER_FIELD_NONZERO;
            break;
        }

        /* (GR_UNABLE: a zero divisor was met, and recorded) */
        if (status == GR_UNABLE && gr_tower_refine(T))
            continue;

        res = GR_TOWER_UNKNOWN;
        break;
    }

    GR_TMP_CLEAR(t, top);
    return res;
}

truth_t
gr_tower_is_zero_at(gr_srcptr x, slong k, gr_tower_t T)
{
    int res = _gr_tower_field_is_zero_at(x, k, T);

    if (res == GR_TOWER_ZERO)
        return T_TRUE;
    if (res == GR_TOWER_NONZERO)
        return T_FALSE;
    if (res == GR_TOWER_UNKNOWN)
        return T_UNKNOWN;

    /* x is a nonzero element of the field; this means that the number
       is nonzero if the transcendental generators are algebraically
       independent, which is only conjectural in general */
    if (!_gr_tower_has_conjectural(T))
        return T_FALSE;

    return _gr_tower_certify_nonzero_at(x, k, T);
}

truth_t
gr_tower_is_zero(gr_srcptr x, gr_tower_t T)
{
    return gr_tower_is_zero_at(x, T->length, T);
}

truth_t
gr_tower_equal_at(gr_srcptr x, gr_srcptr y, slong k, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field_at(T, k);
    gr_ptr t;
    truth_t res;

    GR_TMP_INIT(t, top);

    if (gr_sub(t, x, y, top) == GR_SUCCESS)
        res = gr_tower_is_zero_at(t, k, T);
    else
        res = T_UNKNOWN;

    GR_TMP_CLEAR(t, top);
    return res;
}

truth_t
gr_tower_equal(gr_srcptr x, gr_srcptr y, gr_tower_t T)
{
    return gr_tower_equal_at(x, y, T->length, T);
}

int
gr_tower_inv_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field_at(T, k);
    int status;

    for (;;)
    {
        status = gr_inv(res, x, top);

        /* (GR_UNABLE: a zero divisor was met, and recorded) */
        if (status == GR_UNABLE && gr_tower_refine(T))
            continue;

        return status;
    }
}

int
gr_tower_inv(gr_ptr res, gr_srcptr x, gr_tower_t T)
{
    return gr_tower_inv_at(res, x, T->length, T);
}

int
gr_tower_div_at(gr_ptr res, gr_srcptr x, gr_srcptr y, slong k, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field_at(T, k);
    gr_ptr t;
    int status;

    GR_TMP_INIT(t, top);
    status = gr_tower_inv_at(t, y, k, T);
    if (status == GR_SUCCESS)
        status = gr_mul(res, x, t, top);
    GR_TMP_CLEAR(t, top);
    return status;
}

int
gr_tower_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_tower_t T)
{
    return gr_tower_div_at(res, x, y, T->length, T);
}
