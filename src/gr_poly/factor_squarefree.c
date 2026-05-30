/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011, 2020, 2025, 2026 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* Extract the p-th root of a polynomial over a finite field of characteristic p.
   Precondition: f' = 0. */
static int
_gr_poly_pth_root(gr_poly_t h, const gr_poly_t f, ulong p, gr_ctx_t ctx)
{
    slong deg, i;
    int status = GR_SUCCESS;
    gr_ptr x;

    if (f->length == 0)
        return gr_poly_zero(h, ctx);

    deg = f->length - 1;
    GR_TMP_INIT(x, ctx);

    gr_poly_fit_length(h, deg / p + 1, ctx);

    for (i = 0; i <= (slong)(deg / p); i++)
    {
        status |= gr_poly_get_coeff_scalar(x, f, (slong)(i * p), ctx);
        status |= gr_fq_pth_root(x, x, ctx);
        status |= gr_poly_set_coeff_scalar(h, i, x, ctx);
    }

    _gr_poly_set_length(h, deg / p + 1, ctx);
    _gr_poly_normalise(h, ctx);
    GR_TMP_CLEAR(x, ctx);

    return status;
}

/* Helper before degree-dependent branches */
#define CHECK_PROPER(pol) \
    if (status != GR_SUCCESS || ((pol)->length > 0 && \
        gr_is_zero(gr_poly_coeff_srcptr((pol), \
        (pol)->length - 1, ctx), ctx) != T_FALSE)) \
    { \
        status |= GR_UNABLE; \
        goto cleanup; \
    }

static int
gr_poly_factor_squarefree_finite_field(gr_ptr c, gr_vec_t fac, gr_vec_t exp,
    const gr_poly_t f, gr_ctx_t poly_ctx, gr_ctx_t fmpz_ctx, gr_ctx_t ctx)
{
    gr_poly_t d, t1, v, w, s;
    fmpz_t e, p;
    slong i;
    int status = GR_SUCCESS;

    fmpz_init(e);
    fmpz_init(p);
    gr_poly_init(d, ctx);
    gr_poly_init(t1, ctx);
    gr_poly_init(v, ctx);
    gr_poly_init(w, ctx);
    gr_poly_init(s, ctx);

    status |= gr_poly_get_coeff_scalar(c, f, f->length - 1, ctx);
    status |= gr_poly_derivative(t1, f, ctx);
    CHECK_PROPER(t1)

    /* Case 1: f' = 0 means f = h(x^p); extract h and recurse. */
    if (t1->length == 0)
    {
        gr_vec_t sub_fac, sub_exp;
        gr_ptr sub_c;
        ulong p_ui;
        slong j;

        status |= gr_ctx_fq_prime(p, ctx);
        if (status != GR_SUCCESS)
            goto cleanup;
        p_ui = fmpz_get_ui(p);

        gr_vec_init(sub_fac, 0, poly_ctx);
        gr_vec_init(sub_exp, 0, fmpz_ctx);
        sub_c = gr_heap_init(ctx);

        status |= _gr_poly_pth_root(t1, f, p_ui, ctx);
        status |= gr_poly_factor_squarefree_finite_field(sub_c, sub_fac, sub_exp, t1, poly_ctx, fmpz_ctx, ctx);

        for (j = 0; j < sub_fac->length; j++)
        {
            status |= gr_vec_append(fac,
                gr_vec_entry_ptr(sub_fac, j, poly_ctx), poly_ctx);
            fmpz_mul_ui(e,
                (fmpz *) gr_vec_entry_ptr(sub_exp, j, fmpz_ctx), p_ui);
            status |= gr_vec_append(exp, e, fmpz_ctx);
        }

        gr_vec_clear(sub_fac, poly_ctx);
        gr_vec_clear(sub_exp, fmpz_ctx);
        gr_heap_clear(sub_c, ctx);

        goto cleanup;
    }

    /* Case 2: f' != 0.  Yun's (g_1, g) loop; track the inseparable residual. */
    status |= gr_poly_gcd(d, f, t1, ctx);
    CHECK_PROPER(d)

    if (d->length == 1)
    {
        /* f is already squarefree */
        status |= gr_vec_append(fac, f, poly_ctx);
        fmpz_one(e);
        status |= gr_vec_append(exp, e, fmpz_ctx);
    }
    else
    {
        /* g_1 = f/d,  g = d  (refined at each step) */
        status |= gr_poly_divexact(v, f, d, ctx);
        CHECK_PROPER(v)

        i = 1;
        while (v->length > 1)
        {
            status |= gr_poly_gcd(s, v, d, ctx);
            status |= gr_poly_divexact(w, v, s, ctx);
            CHECK_PROPER(w)

            if (w->length > 1)
            {
                status |= gr_poly_make_monic(w, w, ctx);
                status |= gr_vec_append(fac, w, poly_ctx);
                fmpz_set_ui(e, i);
                status |= gr_vec_append(exp, e, fmpz_ctx);
            }

            status |= gr_poly_divexact(d, d, s, ctx);
            status |= gr_poly_set(v, s, ctx);
            i++;
        }

        CHECK_PROPER(d)

        /* d holds the inseparable part; recurse on d^(1/p). */
        if (d->length > 1)
        {
            gr_vec_t sub_fac, sub_exp;
            gr_ptr sub_c;
            ulong p_ui;
            slong j;

            status |= gr_ctx_fq_prime(p, ctx);
            if (status != GR_SUCCESS)
                goto cleanup;
            p_ui = fmpz_get_ui(p);

            gr_vec_init(sub_fac, 0, poly_ctx);
            gr_vec_init(sub_exp, 0, fmpz_ctx);
            sub_c = gr_heap_init(ctx);

            status |= gr_poly_make_monic(d, d, ctx);
            status |= _gr_poly_pth_root(t1, d, p_ui, ctx);
            status |= gr_poly_factor_squarefree_finite_field(sub_c, sub_fac, sub_exp, t1, poly_ctx, fmpz_ctx, ctx);

            for (j = 0; j < sub_fac->length; j++)
            {
                status |= gr_vec_append(fac,
                    gr_vec_entry_ptr(sub_fac, j, poly_ctx), poly_ctx);
                fmpz_mul_ui(e,
                    (fmpz *) gr_vec_entry_ptr(sub_exp, j, fmpz_ctx), p_ui);
                status |= gr_vec_append(exp, e, fmpz_ctx);
            }

            gr_vec_clear(sub_fac, poly_ctx);
            gr_vec_clear(sub_exp, fmpz_ctx);
            gr_heap_clear(sub_c, ctx);
        }
    }

    for (i = 0; i < fac->length; i++)
        status |= gr_poly_make_monic(gr_vec_entry_ptr(fac, i, poly_ctx),
            gr_vec_entry_srcptr(fac, i, poly_ctx), ctx);

cleanup:
    gr_poly_clear(d, ctx);
    gr_poly_clear(t1, ctx);
    gr_poly_clear(v, ctx);
    gr_poly_clear(w, ctx);
    gr_poly_clear(s, ctx);
    fmpz_clear(e);
    fmpz_clear(p);

    return status;
}

static int
gr_poly_factor_squarefree_ufd_char_0(gr_ptr c, gr_vec_t fac, gr_vec_t exp,
    const gr_poly_t F, gr_ctx_t poly_ctx, gr_ctx_t fmpz_ctx, truth_t is_field, gr_ctx_t ctx)
{
    gr_poly_t f, d, t1, v, w, s;
    gr_ptr lc, lc_pow;
    fmpz_t e;
    slong i, k;
    int status = GR_SUCCESS;

    fmpz_init(e);
    gr_poly_init(f, ctx);
    gr_poly_init(d, ctx);
    gr_poly_init(t1, ctx);
    gr_poly_init(v, ctx);
    gr_poly_init(w, ctx);
    gr_poly_init(s, ctx);
    lc = gr_heap_init(ctx);
    lc_pow = gr_heap_init(ctx);

    if (is_field == T_TRUE)
    {
        status |= gr_poly_get_coeff_scalar(c, F, F->length - 1, ctx);
        status |= gr_poly_make_monic(f, F, ctx);
    }
    else
    {
        /* todo: public vec_content, poly_contents methods */
        status |= gr_poly_get_coeff_scalar(c, F, 0, ctx);
        for (k = 1; k < F->length; k++)
            status |= gr_gcd(c, c, gr_poly_coeff_srcptr(F, k, ctx), ctx);
        status |= gr_poly_divexact_scalar(f, F, c, ctx);
    }

    status |= gr_poly_derivative(t1, f, ctx);
    status |= gr_poly_gcd(d, f, t1, ctx);
    CHECK_PROPER(d)

    if (d->length == 1)
    {
        status |= gr_vec_append(fac, f, poly_ctx);
        fmpz_one(e);
        status |= gr_vec_append(exp, e, fmpz_ctx);
    }
    else
    {
        status |= gr_poly_divexact(v, f, d, ctx);
        status |= gr_poly_divexact(w, t1, d, ctx);

        for (i = 1; ; i++)
        {
            status |= gr_poly_derivative(t1, v, ctx);
            status |= gr_poly_sub(s, w, t1, ctx);
            CHECK_PROPER(s)

            if (s->length == 0)
            {
                CHECK_PROPER(v)

                if (v->length > 1)
                {
                    status |= gr_vec_append(fac, v, poly_ctx);
                    fmpz_set_ui(e, i);
                    status |= gr_vec_append(exp, e, fmpz_ctx);
                }

                break;
            }

            status |= gr_poly_gcd(d, v, s, ctx);
            status |= gr_poly_divexact(v, v, d, ctx);
            status |= gr_poly_divexact(w, s, d, ctx);

            CHECK_PROPER(d)

            if (d->length > 1)
            {
                status |= gr_vec_append(fac, d, poly_ctx);
                fmpz_set_ui(e, i);
                status |= gr_vec_append(exp, e, fmpz_ctx);
            }
        }
    }

    /* Normalize leading coefficients of all factors and adjust unit. */
    if (status == GR_SUCCESS)
    {
        if (is_field == T_TRUE)
        {
            for (i = 0; i < fac->length; i++)
                status |= gr_poly_make_monic(gr_vec_entry_ptr(fac, i, poly_ctx),
                    gr_vec_entry_srcptr(fac, i, poly_ctx), ctx);
        }
        else
        {
            fmpz *expc = exp->entries;
            slong j;

            status |= gr_one(lc, ctx);

            for (j = 0; j < fac->length; j++)
            {
                gr_poly_struct *fac_j = (gr_poly_struct *) gr_vec_entry_ptr(fac, j, poly_ctx);

                status |= gr_poly_canonical_associate(fac_j, NULL, fac_j, ctx);
                CHECK_PROPER(fac_j)
                gr_srcptr lc_j = gr_poly_coeff_srcptr(fac_j, fac_j->length - 1, ctx);
                status |= gr_pow_ui(lc_pow, lc_j, (ulong) expc[j], ctx);
                status |= gr_mul(lc, lc, lc_pow, ctx);
            }

            status |= gr_divexact(lc, gr_poly_coeff_srcptr(f, f->length - 1, ctx), lc, ctx);
            status |= gr_mul(c, c, lc, ctx);
        }
    }

cleanup:
    gr_poly_clear(f, ctx);
    gr_poly_clear(d, ctx);
    gr_poly_clear(t1, ctx);
    gr_poly_clear(v, ctx);
    gr_poly_clear(w, ctx);
    gr_poly_clear(s, ctx);
    gr_heap_clear(lc, ctx);
    gr_heap_clear(lc_pow, ctx);
    fmpz_clear(e);

    return status;
}


int
gr_poly_factor_squarefree(gr_ptr c, gr_vec_t fac, gr_vec_t exp,
    const gr_poly_t F, gr_ctx_t ctx)
{
    gr_ctx_t poly_ctx, fmpz_ctx;
    int status = GR_SUCCESS;
    truth_t is_field, is_fin_char;

    gr_ctx_init_gr_poly(poly_ctx, ctx);
    gr_ctx_init_fmpz(fmpz_ctx);

    gr_vec_set_length(fac, 0, poly_ctx);
    gr_vec_set_length(exp, 0, fmpz_ctx);

    if (F->length == 0)
    {
        status |= gr_zero(c, ctx);
        goto done;
    }

    if (gr_is_zero(gr_poly_coeff_srcptr(F, F->length - 1, ctx), ctx) != T_FALSE)
    {
        status = GR_UNABLE;
        goto done;
    }

    if (F->length == 1)
    {
        status |= gr_set(c, F->coeffs, ctx);
        goto done;
    }

    is_field = gr_ctx_is_field(ctx);
    is_fin_char = gr_ctx_is_finite_characteristic(ctx);

    if (is_field == T_TRUE && is_fin_char == T_TRUE)
    {
        status |= gr_poly_factor_squarefree_finite_field(c, fac, exp, F, poly_ctx, fmpz_ctx, ctx);
    }
    else if (gr_ctx_is_unique_factorization_domain(ctx) == T_TRUE && is_fin_char == T_FALSE)
    {
        status |= gr_poly_factor_squarefree_ufd_char_0(c, fac, exp, F, poly_ctx, fmpz_ctx, is_field, ctx);
    }
    else
    {
        status = GR_UNABLE;
    }

done:
    gr_ctx_clear(poly_ctx);
    gr_ctx_clear(fmpz_ctx);

    return status;
}

