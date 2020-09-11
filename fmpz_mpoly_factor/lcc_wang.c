/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "n_poly.h"


int fmpz_mpoly_factor_lcc_wang(
    fmpz_mpoly_struct * lc_divs,
    const fmpz_mpoly_factor_t lcAfac,
    const fmpz_t Auc,
    const fmpz_poly_struct * Auf,
    slong r,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    fmpz * lcAfaceval = _fmpz_vec_init(lcAfac->num);
    fmpz * d = _fmpz_vec_init(1 + lcAfac->num);
    fmpz * dtilde = _fmpz_vec_init(r);
    fmpz_t Q, R;
    fmpz_mpoly_t t;
    slong N, * offsets, * shifts, * starts, * ends, * stops;
    ulong mask, * es;
    fmpz * T;

    fmpz_init(Q);
    fmpz_init(R);

    fmpz_mpoly_init(t, ctx);

    starts = FLINT_ARRAY_ALLOC(n + 1, slong);
    ends   = FLINT_ARRAY_ALLOC(n + 1, slong);
    stops  = FLINT_ARRAY_ALLOC(n + 1, slong);
    es     = FLINT_ARRAY_ALLOC(n + 1, ulong);
    T      = FLINT_ARRAY_ALLOC(n + 2, fmpz);
    for (i = 0; i < n + 2; i++)
        fmpz_init(T + i);

    offsets = FLINT_ARRAY_ALLOC(n + 1, slong);
    shifts  = FLINT_ARRAY_ALLOC(n + 1, slong);

    for (j = 0; j < lcAfac->num; j++)
    {
        fmpz_mpoly_struct * P = lcAfac->poly + j;

        for (i = 0; i < n + 1; i++)
            mpoly_gen_offset_shift_sp(offsets + i, shifts + i, i, P->bits, ctx->minfo);

        mask = (-UWORD(1)) >> (FLINT_BITS - P->bits);
        N = mpoly_words_per_exp_sp(P->bits, ctx->minfo);
        _fmpz_mpoly_evaluate_rest_fmpz(T, starts, ends, stops, es,
                                      P->coeffs, P->exps, P->length, 1, alpha,
                                              offsets, shifts, N, mask, n + 1);
        fmpz_set(lcAfaceval + j, T + 0);
    }

    fmpz_mul(d + 0, Auc, lcAfac->constant);
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_abs(Q, lcAfaceval + i);
        if (fmpz_is_one(Q) || fmpz_is_zero(Q))
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            fmpz_set(R, d + j);
            while (!fmpz_is_one(R))
            {
                fmpz_gcd(R, R, Q);
                fmpz_divexact(Q, Q, R);
                if (fmpz_is_one(Q))
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        fmpz_set(d + i + 1, Q);
    }

    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_one(lc_divs + j, ctx);
        fmpz_one(dtilde + j);
        fmpz_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            fmpz_abs(Q, lcAfaceval + i);
            if (fmpz_cmp_ui(Q, 2) < 0)
                continue;
            k = fmpz_remove(R, R, Q);
            fmpz_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            fmpz_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
            fmpz_pow_ui(Q, lcAfaceval + i, k);
            fmpz_mul(dtilde + j, dtilde + j, Q);
        }
    }

    for (j = 0; j < r; j++)
    {
        FLINT_ASSERT(Auf[j].length > 0);
        fmpz_gcd(T + 0, Auf[j].coeffs + Auf[j].length - 1, dtilde + j);
        fmpz_fdiv_qr(Q, R, Auf[j].coeffs + Auf[j].length - 1, T + 0);
        if (!fmpz_is_zero(R))
        {
            success = 0;
            goto cleanup;
        }
        fmpz_mpoly_scalar_mul_fmpz(lc_divs + j, lc_divs + j, Q, ctx);
    }

    success = 1;

cleanup:

    fmpz_clear(Q);
    fmpz_clear(R);
    fmpz_mpoly_clear(t, ctx);
    _fmpz_vec_clear(lcAfaceval, lcAfac->num);
    _fmpz_vec_clear(d, 1 + lcAfac->num);
    _fmpz_vec_clear(dtilde, r);

    for (i = 0; i < n + 2; i++)
        fmpz_clear(T + i);
    flint_free(T);
    flint_free(starts);
    flint_free(ends);
    flint_free(stops);
    flint_free(es);

    flint_free(offsets);
    flint_free(shifts);

    return success;
}

