/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "fq_nmod_mpoly_factor.h"


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
    n_poly_struct * d;
    n_poly_t Q, R;
    fq_nmod_mpoly_t t;
    slong N, * offsets, * shifts, * starts, * ends, * stops;
    ulong mask, * es;
    n_poly_struct * T;

    n_poly_init(Q);
    n_poly_init(R);
    fq_nmod_mpoly_init(t, ctx);

    lcAfaceval = FLINT_ARRAY_ALLOC(lcAfac->num, n_poly_struct);
    for (i = 0; i < lcAfac->num; i++)
        n_poly_init(lcAfaceval + i);

    d = FLINT_ARRAY_ALLOC(lcAfac->num + 1, n_poly_struct);
    for (i = 0; i < lcAfac->num + 1; i++)
        n_poly_init(d + i);

    starts = FLINT_ARRAY_ALLOC(n + 1, slong);
    ends   = FLINT_ARRAY_ALLOC(n + 1, slong);
    stops  = FLINT_ARRAY_ALLOC(n + 1, slong);
    es     = FLINT_ARRAY_ALLOC(n + 1, ulong);
    T      = FLINT_ARRAY_ALLOC(n + 2, n_poly_struct);
    for (i = 0; i < n + 2; i++)
        n_poly_init(T + i);

    offsets = FLINT_ARRAY_ALLOC(n + 1, slong);
    shifts  = FLINT_ARRAY_ALLOC(n + 1, slong);

    /* init done */

    for (j = 0; j < lcAfac->num; j++)
    {
        fq_nmod_mpoly_struct * P = lcAfac->poly + j;

        for (i = 0; i < n + 1; i++)
            mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, P->bits, ctx->minfo);

        mask = (-UWORD(1)) >> (FLINT_BITS - P->bits);
        N = mpoly_words_per_exp_sp(P->bits, ctx->minfo);
        _fq_nmod_mpoly_eval_rest_n_fq_poly(T, starts, ends, stops, es,
                                      P->coeffs, P->exps, P->length, 1, alpha,
                                  offsets, shifts, N, mask, n + 1, ctx->fqctx);

        n_fq_poly_set(lcAfaceval + j, T + 0, ctx->fqctx);
    }

    n_fq_poly_set(d + 0, Auc, ctx->fqctx);
    for (i = 0; i < lcAfac->num; i++)
    {
        n_fq_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
        if (n_poly_degree(Q) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            n_fq_poly_set(R, d + j, ctx->fqctx);
            while (n_poly_degree(R) > 0)
            {
                n_fq_poly_gcd(R, R, Q, ctx->fqctx);
                n_fq_poly_divrem(Q, T + 0, Q, R, ctx->fqctx);
                FLINT_ASSERT(n_poly_is_zero(T + 0));
                if (n_poly_degree(Q) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        n_fq_poly_set(d + i + 1, Q, ctx->fqctx);
    }

    for (j = 0; j < r; j++)
    {
        fq_nmod_mpoly_one(lc_divs + j, ctx);
        n_fq_poly_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->fqctx);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            n_fq_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
            if (n_poly_degree(Q) < 1)
                continue;
            k = n_fq_poly_remove(R, Q, ctx->fqctx);
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

    for (i = 0; i < n + 2; i++)
        n_poly_clear(T + i);
    flint_free(T);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);

    return success;
}
