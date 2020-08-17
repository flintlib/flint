/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


int fq_zech_mpoly_factor_lcc_wang(
    fq_zech_mpoly_struct * lc_divs,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_poly_t Auc,
    const fq_zech_bpoly_struct * Auf,
    slong r,
    const fq_zech_poly_struct * alpha,
    const fq_zech_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fq_zech_poly_struct * lcAfaceval;
    fq_zech_poly_struct * d;
    fq_zech_poly_t Q, R;
    fq_zech_mpoly_t t;
    slong N, * offsets, * shifts, * starts, * ends, * stops;
    ulong mask, * es;
    fq_zech_poly_struct * T;

    fq_zech_poly_init(Q, ctx->fqctx);
    fq_zech_poly_init(R, ctx->fqctx);
    fq_zech_mpoly_init(t, ctx);

    lcAfaceval = FLINT_ARRAY_ALLOC(lcAfac->num, fq_zech_poly_struct);
    for (i = 0; i < lcAfac->num; i++)
        fq_zech_poly_init(lcAfaceval + i, ctx->fqctx);

    d = FLINT_ARRAY_ALLOC(lcAfac->num + 1, fq_zech_poly_struct);
    for (i = 0; i < lcAfac->num + 1; i++)
        fq_zech_poly_init(d + i, ctx->fqctx);

    starts = FLINT_ARRAY_ALLOC(n + 1, slong);
    ends   = FLINT_ARRAY_ALLOC(n + 1, slong);
    stops  = FLINT_ARRAY_ALLOC(n + 1, slong);
    es     = FLINT_ARRAY_ALLOC(n + 1, ulong);
    T      = FLINT_ARRAY_ALLOC(n + 2, fq_zech_poly_struct);
    for (i = 0; i < n + 2; i++)
        fq_zech_poly_init(T + i, ctx->fqctx);

    offsets = FLINT_ARRAY_ALLOC(n + 1, slong);
    shifts  = FLINT_ARRAY_ALLOC(n + 1, slong);

    /* init done */

    for (j = 0; j < lcAfac->num; j++)
    {
        fq_zech_mpoly_struct * P = lcAfac->poly + j;

        for (i = 0; i < n + 1; i++)
            mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, P->bits, ctx->minfo);

        mask = (-UWORD(1)) >> (FLINT_BITS - P->bits);
        N = mpoly_words_per_exp_sp(P->bits, ctx->minfo);
        _fq_zech_mpoly_eval_rest_fq_zech_poly(T, starts, ends, stops, es,
                                      P->coeffs, P->exps, P->length, 1, alpha,
                                  offsets, shifts, N, mask, n + 1, ctx->fqctx);

        fq_zech_poly_set(lcAfaceval + j, T + 0, ctx->fqctx);
    }

    fq_zech_poly_set(d + 0, Auc, ctx->fqctx);
    for (i = 0; i < lcAfac->num; i++)
    {
        fq_zech_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
        if (fq_zech_poly_degree(Q, ctx->fqctx) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fq_zech_poly_set(R, d + j, ctx->fqctx);
            while (fq_zech_poly_degree(R, ctx->fqctx) > 0)
            {
                fq_zech_poly_gcd(R, R, Q, ctx->fqctx);
                fq_zech_poly_divrem(Q, T + 0, Q, R, ctx->fqctx);
                FLINT_ASSERT(fq_zech_poly_is_zero(T + 0, ctx->fqctx));
                if (fq_zech_poly_degree(Q, ctx->fqctx) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fq_zech_poly_set(d + i + 1, Q, ctx->fqctx);
    }

    for (j = 0; j < r; j++)
    {
        fq_zech_mpoly_one(lc_divs + j, ctx);
        fq_zech_poly_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->fqctx);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fq_zech_poly_make_monic(Q, lcAfaceval + i, ctx->fqctx);
            if (fq_zech_poly_degree(Q, ctx->fqctx) < 1)
                continue;
            k = fq_zech_poly_remove(R, Q, ctx->fqctx);
            fq_zech_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fq_zech_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    fq_zech_poly_clear(Q, ctx->fqctx);
    fq_zech_poly_clear(R, ctx->fqctx);
    fq_zech_mpoly_clear(t, ctx);

    for (i = 0; i < lcAfac->num; i++)
        fq_zech_poly_clear(lcAfaceval + i, ctx->fqctx);
    flint_free(lcAfaceval);

    for (i = 0; i < lcAfac->num + 1; i++)
        fq_zech_poly_clear(d + i, ctx->fqctx);
    flint_free(d);

    for (i = 0; i < n + 2; i++)
        fq_zech_poly_clear(T + i, ctx->fqctx);
    flint_free(T);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);

    return success;
}
