/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


ulong n_poly_fq_remove(
    n_poly_t f,
    const n_poly_t g,
    const fq_nmod_ctx_t ctx)
{
    n_poly_t q, r;
    ulong i = 0;

    n_poly_init(q);
    n_poly_init(r);

    while (1)
    {
        if (f->length < g->length)
            break;
        n_poly_fq_divrem(q, r, f, g, ctx);
        if (r->length == 0)
            n_poly_swap(q, f);
        else
            break;
        i++;
    }

    n_poly_clear(q);
    n_poly_clear(r);

    return i;
}

int fq_nmod_mpoly_evaluate_all_n_poly_fq(
    n_poly_t A,
    const fq_nmod_mpoly_t B,
    const n_poly_struct * C,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fq_nmod_poly_t a;
    fq_nmod_poly_struct * t1, ** t2;

    t1 = FLINT_ARRAY_ALLOC(nvars, fq_nmod_poly_struct);
    t2 = FLINT_ARRAY_ALLOC(nvars, fq_nmod_poly_struct *);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_poly_init(t1 + i, ctx->fqctx);
        n_poly_fq_get_fq_nmod_poly(t1 + i, C + i, ctx->fqctx);
        t2[i] = t1 + i;
    }

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_mpoly_compose_fq_nmod_poly(a, B, t2, ctx);
    n_poly_fq_set_fq_nmod_poly(A, a, ctx->fqctx);
    fq_nmod_poly_clear(a, ctx->fqctx);

    for (i = 0; i < nvars; i++)
        fq_nmod_poly_clear(t1 + i, ctx->fqctx);

    flint_free(t1);
    flint_free(t2);

    return 1;
}


int fq_nmod_mpoly_factor_lcc_wang(
    fq_nmod_mpoly_struct * lc_divs,
    const fq_nmod_mpoly_factor_t lcAfac,
    const n_poly_t Auc,
    const n_bpoly_struct * Auf,
    slong r,
    const n_poly_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    n_poly_struct * lcAfaceval;
    n_poly_struct * salpha;
    n_poly_struct * d;
    n_poly_t Q, R;
    fq_nmod_mpoly_t t;

    n_poly_init(Q);
    n_poly_init(R);
    fq_nmod_mpoly_init(t, ctx);

    salpha = FLINT_ARRAY_ALLOC((n + 1), n_poly_struct);
    n_poly_init(salpha + 0);
    for (i = 0; i < n; i++)
        salpha[i + 1] = alpha[i];

    lcAfaceval = FLINT_ARRAY_ALLOC(lcAfac->num, n_poly_struct);
    for (i = 0; i < lcAfac->num; i++)
        n_poly_init(lcAfaceval + i);

    d = FLINT_ARRAY_ALLOC(lcAfac->num + 1, n_poly_struct);
    for (i = 0; i < lcAfac->num + 1; i++)
        n_poly_init(d + i);

    /* init done */

    for (j = 0; j < lcAfac->num; j++)
        fq_nmod_mpoly_evaluate_all_n_poly_fq(lcAfaceval + j, lcAfac->poly + j, salpha, ctx);

    n_poly_fq_set(d + 0, Auc, ctx->fqctx);
    for (i = 0; i < lcAfac->num; i++)
    {
        n_poly_fq_make_monic(Q, lcAfaceval + i, ctx->fqctx);
        if (n_poly_degree(Q) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            n_poly_fq_set(R, d + j, ctx->fqctx);
            while (n_poly_degree(R) > 0)
            {
                n_poly_fq_gcd(R, R, Q, ctx->fqctx);
                n_poly_fq_divrem(Q, salpha + 0, Q, R, ctx->fqctx);
                FLINT_ASSERT(n_poly_is_zero(salpha + 0));
                if (n_poly_degree(Q) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        n_poly_fq_set(d + i + 1, Q, ctx->fqctx);
    }

    for (j = 0; j < r; j++)
    {
        fq_nmod_mpoly_one(lc_divs + j, ctx);
        n_poly_fq_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->fqctx);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            n_poly_fq_make_monic(Q, lcAfaceval + i, ctx->fqctx);
            if (n_poly_degree(Q) < 1)
                continue;
            k = n_poly_fq_remove(R, Q, ctx->fqctx);
            fq_nmod_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fq_nmod_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    n_poly_clear(Q);
    n_poly_clear(R);
    fq_nmod_mpoly_clear(t, ctx);

    for (i = 0; i < lcAfac->num; i++)
        n_poly_clear(lcAfaceval + i);
    flint_free(lcAfaceval);

    for (i = 0; i < lcAfac->num + 1; i++)
        n_poly_clear(d + i);
    flint_free(d);

    n_poly_clear(salpha + 0);
    flint_free(salpha);

    return success;
}
