/*
    Copyright (C) 2018,2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           ulong k, const nmod_mpoly_ctx_t ctx)
{
    slong i, exp_bits, N;
    fmpz * maxBfields;
    ulong * cmpmask;
    ulong * Bexps;
    int freeBexps;
    nmod_mpoly_t T, Atemp;
    nmod_mpoly_struct * R;

    TMP_INIT;

    if (k == 0)
    {
        nmod_mpoly_set_ui(A, ctx->mod.n > 1, ctx);
        return 1;
    }

    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return 1;
    }

    if (k == 1)
    {
        nmod_mpoly_set(A, B, ctx);
        return 1;
    }

    if (k == 2)
    {
        nmod_mpoly_mul(A, B, B, ctx);
        return 1;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);

    _fmpz_vec_scalar_mul_ui(maxBfields, maxBfields, ctx->minfo->nfields, k);

    exp_bits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);

    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    N = mpoly_words_per_exp(exp_bits, ctx->minfo);

    if (B->length == 1)
    {
        /* powering a monomial */
        nmod_mpoly_fit_length_reset_bits(A, 1, exp_bits, ctx);

        if (B->bits == exp_bits && B != A)
            mpoly_monomial_mul_ui_mp(A->exps, B->exps, N, k);
        else
            mpoly_pack_vec_fmpz(A->exps, maxBfields, exp_bits,
                                                       ctx->minfo->nfields, 1);

        A->coeffs[0] = nmod_pow_ui(B->coeffs[0], k, ctx->mod);

        _nmod_mpoly_set_length(A, A->coeffs[0] != 0, ctx);

        goto cleanup;
    }

    freeBexps = 0;
    Bexps = B->exps;
    if (exp_bits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, exp_bits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    if (A == B)
    {
        nmod_mpoly_init3(Atemp, B->length, exp_bits, ctx);
        R = Atemp;
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(A, B->length, exp_bits, ctx);
        R = A;
    }

    nmod_mpoly_init3(T, B->length, exp_bits, ctx);

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    if (ctx->mod.n > 99999 || !n_is_prime(ctx->mod.n))
    {
        _nmod_mpoly_pow_rmul(R, B->coeffs, Bexps, B->length, k,
                                                      N, cmpmask, ctx->mod, T);
    }
    else
    {
        ulong ne, kmodn;
        nmod_mpoly_t S;

        nmod_mpoly_init3(S, B->length, exp_bits, ctx);

        mpoly_monomial_zero(R->exps, N);
        R->coeffs[0] = 1;
        R->length = 1;

        for (ne = 1; k > 0; k = k/ctx->mod.n, ne = ne * ctx->mod.n)
        {
            NMOD_RED(kmodn, k, ctx->mod);

            /* R *= B^(n^e*(k%n)) */

            if (kmodn == 0)
                continue;

            _nmod_mpoly_pow_rmul(S, B->coeffs, Bexps, B->length, kmodn,
                                                      N, cmpmask, ctx->mod, T);

            mpoly_monomial_mul_ui_mp(S->exps, S->exps, N*S->length, ne);

            if (nmod_mpoly_is_one(R, ctx))
            {
                nmod_mpoly_swap(R, S, ctx);
            }
            else
            {
                _nmod_mpoly_mul_johnson(T, R->coeffs, R->exps, R->length,
                                           S->coeffs, S->exps, S->length,
                                           exp_bits, N, cmpmask, ctx->mod);
                nmod_mpoly_swap(R, T, ctx);
            }
        }

        nmod_mpoly_clear(S, ctx);
    }

    nmod_mpoly_clear(T, ctx);

    if (A == B)
    {
        nmod_mpoly_swap(A, Atemp, ctx);
        nmod_mpoly_clear(Atemp, ctx);
    }

    if (freeBexps)
        flint_free(Bexps);

cleanup:

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;

    return 1;
}

