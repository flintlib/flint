/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"


void nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           slong k, const nmod_mpoly_ctx_t ctx)
{
    slong i, exp_bits, N;
    fmpz * Bmaxfields;
    ulong * Bexps;
    ulong * cmpmask;
    int freeBexps;
    TMP_INIT;

    if (k == 0)
    {
        nmod_mpoly_set_ui(A, 1, ctx);
        return;
    }

    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    if (k == 1)
    {
        nmod_mpoly_set(A, B, ctx);
        return;
    }

    if (k == 2)
    {
        nmod_mpoly_mul_johnson(A, B, B, ctx);
        return;
    }

    TMP_START;

    Bmaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(Bmaxfields + i);
    }
    mpoly_max_fields_fmpz(Bmaxfields, B->exps, B->length, B->bits, ctx->minfo);

    _fmpz_vec_scalar_mul_ui(Bmaxfields, Bmaxfields, ctx->minfo->nfields, k);

    exp_bits = _fmpz_vec_max_bits(Bmaxfields, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);

    freeBexps = 0;
    Bexps = B->exps;
    if (exp_bits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, exp_bits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    if (B->length == 1)
    {
        /* powering a monomial */
        nmod_mpoly_fit_length(A, 1, ctx);
        nmod_mpoly_fit_bits(A, exp_bits, ctx);
        A->bits = exp_bits;

        if (exp_bits <= FLINT_BITS)
        {
            mpoly_monomial_mul_si(A->exps, Bexps, N, k);
        } else
        {
            mpoly_monomial_mul_ui_mp(A->exps, Bexps, N, k);
        }

        A->coeffs[0] = nmod_pow_ui(B->coeffs[0], k, ctx->ffinfo->mod);
        _nmod_mpoly_set_length(A, WORD(A->coeffs[0] != 0), ctx);

    } else {

        cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

        if (1 || !n_is_prime(ctx->ffinfo->mod.n))
        {
            fmpz_mpoly_t T;
            slong Alen, Tlen;
            fmpz * Bcoeffs_fmpz;

            T->alloc = (k*(B->length - 1) + 1);
            T->coeffs = (fmpz *) flint_calloc(T->alloc, sizeof(fmpz));
            T->exps   = (ulong *) flint_malloc(T->alloc*N*sizeof(ulong));
            T->length = 0;
            T->bits = exp_bits;

            Bcoeffs_fmpz = _fmpz_vec_init(B->length);
            _fmpz_vec_set_nmod_vec(Bcoeffs_fmpz, B->coeffs, B->length, ctx->ffinfo->mod);

            Tlen = _fmpz_mpoly_pow_fps(&T->coeffs, &T->exps, &T->alloc,
                       Bcoeffs_fmpz, Bexps, B->length, k, exp_bits, N, cmpmask);

            _fmpz_vec_clear(Bcoeffs_fmpz, B->length);

            nmod_mpoly_fit_length(A, Tlen, ctx);
            nmod_mpoly_fit_bits(A, exp_bits, ctx);
            A->bits = exp_bits;

            Alen = 0;
            for (i = 0; i < Tlen; i++)
            {
                A->coeffs[Alen] = fmpz_fdiv_ui(T->coeffs + i, ctx->ffinfo->mod.n);
                mpoly_monomial_set(A->exps + N*Alen, T->exps + N*i, N);
                Alen += (A->coeffs[Alen] != UWORD(0));
            }

            _nmod_mpoly_set_length(A, Alen, ctx);
      
        } else
        {
            assert(0);
/*
            nmod_mpoly_t r;
            nmod_mpoly_one(r, ctx);
            while (k > 0)
            {
                nmod_mpoly_pow_fps(t, a, k % n, ctx);
                nmod_mpoly_mul(r, r, t);
                nmod_mpoly_pth_power(a, ctx);
                k = k / n;
            }
            nmod_mpoly_swap(r, poly1);
*/
        }
    }

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(Bmaxfields + i);
    }

    if (freeBexps)
    {
        flint_free(Bexps);
    }

    TMP_END;
    return;
}
