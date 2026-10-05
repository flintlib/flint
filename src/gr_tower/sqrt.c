/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Square roots inside a tower of quadratic extensions, without lattice
    reduction: with F_k = F_{k-1}(a), a^2 = r (after completing the
    square), x = u + v a is a square in F_k iff N = u^2 - r v^2 is a
    square s in F_{k-1} and (u + s)/2 or (u - s)/2 is a square p^2 in
    F_{k-1}; then sqrt(x) = p + v/(2p) a. This recurses down to the base
    field. Steps of degree other than 2 are not handled (GR_UNABLE).
*/

#include "fmpq.h"
#include "fmpz_mpoly_q.h"
#include "acb.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/* Square root in the base field: a rational number, or a rational
   function of the transcendental generators with square numerator and
   denominator. */
static int
_base_sqrt(gr_ptr res, gr_srcptr x, gr_tower_t T)
{
    gr_ctx_struct * base = T->base;

    if (base->which_ring == GR_CTX_FMPQ)
    {
        const fmpq * q = x;
        fmpq * r = res;

        if (fmpq_sgn(q) < 0)
            return GR_DOMAIN;
        if (!fmpz_is_square(fmpq_numref(q)) || !fmpz_is_square(fmpq_denref(q)))
            return GR_DOMAIN;
        fmpz_sqrt(fmpq_numref(r), fmpq_numref(q));
        fmpz_sqrt(fmpq_denref(r), fmpq_denref(q));
        return GR_SUCCESS;
    }
    else if (base->which_ring == GR_CTX_FMPZ_MPOLY_Q)
    {
        const fmpz_mpoly_q_struct * q = x;
        fmpz_mpoly_q_struct * r = res;
        const fmpz_mpoly_ctx_struct * mctx = gr_ctx_fmpz_mpoly_q_mctx(base);
        fmpz_mpoly_t a, b;
        int ok;

        fmpz_mpoly_init(a, mctx);
        fmpz_mpoly_init(b, mctx);
        ok = fmpz_mpoly_sqrt(a, fmpz_mpoly_q_numref(q), mctx) && fmpz_mpoly_sqrt(b, fmpz_mpoly_q_denref(q), mctx);
        if (ok)
        {
            fmpz_mpoly_swap(fmpz_mpoly_q_numref(r), a, mctx);
            fmpz_mpoly_swap(fmpz_mpoly_q_denref(r), b, mctx);
            fmpz_mpoly_q_canonicalise(r, mctx);
        }
        fmpz_mpoly_clear(a, mctx);
        fmpz_mpoly_clear(b, mctx);
        return ok ? GR_SUCCESS : GR_DOMAIN;
    }

    return GR_UNABLE;
}

/*
    Some square root of x in F_k (sign unspecified). The recursion through
    the quadratic steps branches (up to three square roots one level
    down), so its cost grows exponentially with the number of quadratic
    steps above the level of x; *budget bounds the number of calls (then
    GR_UNABLE: the caller adjoins a root instead; GR_TOWER_OPT_SQRT_BUDGET).
*/
static int
_sqrt_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T, slong * budget)
{
    gr_ctx_struct * Fk, * below;
    gr_ptr u, v, r, N, s, p, q, t, half_b;
    const gr_poly_struct * m;
    int status;

    if (--(*budget) < 0)
        return GR_UNABLE;

    if (k == 0)
        return _base_sqrt(res, x, T);

    if (gr_tower_step_degree(T, k) != 2)
        return GR_UNABLE;

    Fk = gr_tower_field_at(T, k);
    below = gr_tower_field_at(T, k - 1);
    m = gr_tower_step_minpoly(T, k);   /* X^2 + b X + c */

    GR_TMP_INIT5(u, v, r, N, s, below);
    GR_TMP_INIT4(p, q, t, half_b, below);

    /* x = u + v a; with a' = a + b/2, a'^2 = r = b^2/4 - c, x = (u - v b/2) + v a' */
    {
        gr_poly_t xp;
        gr_poly_init(xp, below);
        status = gr_poly_quotient_get_poly(xp, x, Fk);
        if (xp->length >= 1)
            status |= gr_set(u, gr_poly_coeff_srcptr(xp, 0, below), below);
        else
            status |= gr_zero(u, below);
        if (xp->length >= 2)
            status |= gr_set(v, gr_poly_coeff_srcptr(xp, 1, below), below);
        else
            status |= gr_zero(v, below);
        gr_poly_clear(xp, below);
    }

    status |= gr_div_ui(half_b, gr_poly_coeff_srcptr(m, 1, below), 2, below);
    status |= gr_sqr(r, half_b, below);
    status |= gr_sub(r, r, gr_poly_coeff_srcptr(m, 0, below), below);
    status |= gr_mul(t, v, half_b, below);
    status |= gr_sub(u, u, t, below);

    if (status != GR_SUCCESS)
        goto cleanup;

    if (gr_tower_is_zero_at(v, k - 1, T) == T_TRUE)
    {
        /* x = u: sqrt in F_{k-1}, or q a' with q^2 r = u */
        status = _sqrt_at(p, u, k - 1, T, budget);
        if (status == GR_SUCCESS)
        {
            status = gr_set_other(res, p, below, Fk);
            goto cleanup;
        }
        if (status != GR_DOMAIN)
            goto cleanup;

        status = gr_tower_div_at(t, u, r, k - 1, T);
        if (status != GR_SUCCESS)
            goto cleanup;
        status = _sqrt_at(q, t, k - 1, T, budget);
        if (status != GR_SUCCESS)
            goto cleanup;
        status = gr_zero(p, below);
    }
    else
    {
        /* N = u^2 - r v^2 must be a square s */
        status = gr_sqr(N, u, below);
        status |= gr_sqr(t, v, below);
        status |= gr_mul(t, t, r, below);
        status |= gr_sub(N, N, t, below);
        if (status != GR_SUCCESS)
            goto cleanup;

        status = _sqrt_at(s, N, k - 1, T, budget);
        if (status != GR_SUCCESS)
            goto cleanup;

        /* p^2 = (u + s)/2 or (u - s)/2 */
        status = gr_add(t, u, s, below);
        status |= gr_div_ui(t, t, 2, below);
        if (status != GR_SUCCESS)
            goto cleanup;
        status = _sqrt_at(p, t, k - 1, T, budget);
        if (status == GR_DOMAIN || (status == GR_SUCCESS && gr_tower_is_zero_at(p, k - 1, T) == T_TRUE))
        {
            status = gr_sub(t, u, s, below);
            status |= gr_div_ui(t, t, 2, below);
            if (status != GR_SUCCESS)
                goto cleanup;
            status = _sqrt_at(p, t, k - 1, T, budget);
        }
        if (status != GR_SUCCESS)
            goto cleanup;

        if (gr_tower_is_zero_at(p, k - 1, T) == T_TRUE)
        {
            status = GR_DOMAIN;
            goto cleanup;
        }

        /* q = v / (2 p) */
        status = gr_mul_ui(t, p, 2, below);
        status |= gr_tower_div_at(q, v, t, k - 1, T);
        if (status != GR_SUCCESS)
            goto cleanup;
    }

    /* res = p + q a' = (p + q b/2) + q a */
    {
        gr_poly_t rp;
        gr_poly_init(rp, below);
        status = gr_mul(t, q, half_b, below);
        status |= gr_add(t, t, p, below);
        status |= gr_poly_set_coeff_scalar(rp, 0, t, below);
        status |= gr_poly_set_coeff_scalar(rp, 1, q, below);
        status |= gr_poly_quotient_set_poly(res, rp, Fk);
        gr_poly_clear(rp, below);
    }

cleanup:
    GR_TMP_CLEAR5(u, v, r, N, s, below);
    GR_TMP_CLEAR4(p, q, t, half_b, below);
    return status;
}

/*
    The principal square root of x in F_k, if it lies there. Returns
    GR_DOMAIN if x is not a square in F_k, GR_UNABLE if the tower has a
    step of degree other than 2 below k (or the sign cannot be decided).
*/
/* (the principal square root of x, and +/- s, for the sign selection) */
typedef struct { gr_srcptr x; gr_srcptr s; slong k; gr_tower_struct * T; } _sqrt_sign_arg;

static int
_sqrt_sign_get(acb_t z, int which, slong prec, void * arg)
{
    _sqrt_sign_arg * a = (_sqrt_sign_arg *) arg;
    int st;
    if (which < 0)
    {
        st = gr_tower_get_acb_at(z, a->x, a->k, prec, a->T);
        if (st == GR_SUCCESS)
            acb_sqrt(z, z, prec);
    }
    else
    {
        st = gr_tower_get_acb_at(z, a->s, a->k, prec, a->T);
        if (st == GR_SUCCESS && which == 1)
            acb_neg(z, z);
    }
    return st;
}

int
gr_tower_sqrt_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
{
    gr_ctx_struct * Fk = gr_tower_field_at(T, k);
    gr_ptr s;
    int status;
    slong budget = GR_TOWER_OPTION(T, GR_TOWER_OPT_SQRT_BUDGET);

    GR_TMP_INIT(s, Fk);
    status = _sqrt_at(s, x, k, T, &budget);
    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR(s, Fk);
        return status;
    }

    if (gr_tower_is_zero_at(s, k, T) == T_TRUE)
    {
        status = gr_zero(res, Fk);
        GR_TMP_CLEAR(s, Fk);
        return status;
    }

    /* the root with the correct sign: +s or -s, by the enclosures */
    {
        _sqrt_sign_arg a = { x, s, k, T };
        int which = _gr_tower_select_candidate(_sqrt_sign_get, &a, T);
        if (which == 0)
            status = gr_set(res, s, Fk);
        else if (which == 1)
            status = gr_neg(res, s, Fk);
        else
            status = GR_UNABLE;
    }

    GR_TMP_CLEAR(s, Fk);
    return status;
}

int
gr_tower_sqrt(gr_ptr res, gr_srcptr x, gr_tower_t T)
{
    return gr_tower_sqrt_at(res, x, T->length, T);
}
