/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "fmpz_mod_mpoly_factor.h"


int fmpz_mod_mpoly_factor_lcc_wang(
    fmpz_mod_mpoly_struct * lc_divs,
    const fmpz_mod_mpoly_factor_t lcAfac,
    const fmpz_mod_poly_t Auc,
    const fmpz_mod_bpoly_struct * Auf,
    slong r,
    const fmpz_mod_poly_struct * alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fmpz_mod_poly_struct * lcAfaceval;
    fmpz_mod_poly_struct * d;
    fmpz_mod_poly_t Q, R;
    fmpz_mod_mpoly_t t;
    slong N, * offsets, * shifts, * starts, * ends, * stops;
    ulong mask, * es;
    fmpz_mod_poly_struct * T;

    fmpz_mod_poly_init(Q, ctx->ffinfo);
    fmpz_mod_poly_init(R, ctx->ffinfo);
    fmpz_mod_mpoly_init(t, ctx);

    lcAfaceval = FLINT_ARRAY_ALLOC(lcAfac->num, fmpz_mod_poly_struct);
    for (i = 0; i < lcAfac->num; i++)
        fmpz_mod_poly_init(lcAfaceval + i, ctx->ffinfo);

    d = FLINT_ARRAY_ALLOC(lcAfac->num + 1, fmpz_mod_poly_struct);
    for (i = 0; i < lcAfac->num + 1; i++)
        fmpz_mod_poly_init(d + i, ctx->ffinfo);

    starts = FLINT_ARRAY_ALLOC(n + 1, slong);
    ends   = FLINT_ARRAY_ALLOC(n + 1, slong);
    stops  = FLINT_ARRAY_ALLOC(n + 1, slong);
    es     = FLINT_ARRAY_ALLOC(n + 1, ulong);
    T      = FLINT_ARRAY_ALLOC(n + 2, fmpz_mod_poly_struct);
    for (i = 0; i < n + 2; i++)
        fmpz_mod_poly_init(T + i, ctx->ffinfo);

    offsets = FLINT_ARRAY_ALLOC(n + 1, slong);
    shifts  = FLINT_ARRAY_ALLOC(n + 1, slong);

    for (j = 0; j < lcAfac->num; j++)
    {
        fmpz_mod_mpoly_struct * P = lcAfac->poly + j;

        for (i = 0; i < n + 1; i++)
            mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, P->bits, ctx->minfo);

        mask = (-UWORD(1)) >> (FLINT_BITS - P->bits);
        N = mpoly_words_per_exp_sp(P->bits, ctx->minfo);
        _fmpz_mod_mpoly_evaluate_rest_fmpz_mod_poly(T, starts, ends, stops, es,
                                      P->coeffs, P->exps, P->length, 1, alpha,
                                 offsets, shifts, N, mask, n + 1, ctx->ffinfo);

        fmpz_mod_poly_set(lcAfaceval + j, T + 0, ctx->ffinfo);
    }

    fmpz_mod_poly_set(d + 0, Auc, ctx->ffinfo);
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_mod_poly_make_monic(Q, lcAfaceval + i, ctx->ffinfo);
        if (fmpz_mod_poly_degree(Q, ctx->ffinfo) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fmpz_mod_poly_set(R, d + j, ctx->ffinfo);
            while (fmpz_mod_poly_degree(R, ctx->ffinfo) > 0)
            {
                fmpz_mod_poly_gcd(R, R, Q, ctx->ffinfo);
                fmpz_mod_poly_divrem(Q, T + 0, Q, R, ctx->ffinfo);
                if (fmpz_mod_poly_degree(Q, ctx->ffinfo) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fmpz_mod_poly_set(d + i + 1, Q, ctx->ffinfo);
    }

    for (j = 0; j < r; j++)
    {
        fmpz_mod_mpoly_one(lc_divs + j, ctx);
        fmpz_mod_poly_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->ffinfo);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fmpz_mod_poly_make_monic(Q, lcAfaceval + i, ctx->ffinfo);
            if (fmpz_mod_poly_degree(Q, ctx->ffinfo) < 1)
                continue;
            k = fmpz_mod_poly_remove(R, Q, ctx->ffinfo);
            fmpz_mod_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fmpz_mod_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    fmpz_mod_poly_clear(Q, ctx->ffinfo);
    fmpz_mod_poly_clear(R, ctx->ffinfo);
    fmpz_mod_mpoly_clear(t, ctx);

    for (i = 0; i < lcAfac->num; i++)
        fmpz_mod_poly_clear(lcAfaceval + i, ctx->ffinfo);
    flint_free(lcAfaceval);

    for (i = 0; i < lcAfac->num + 1; i++)
        fmpz_mod_poly_clear(d + i, ctx->ffinfo);
    flint_free(d);

    for (i = 0; i < n + 2; i++)
        fmpz_mod_poly_clear(T + i, ctx->ffinfo);
    flint_free(T);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);

    return success;
}
