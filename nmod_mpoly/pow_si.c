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


void nmod_mpoly_pow_si(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           slong k, const nmod_mpoly_ctx_t ctx)
{
    slong i, exp_bits, N;
    fmpz * Bmaxfields;
    ulong * Bexps;
    ulong * cmpmask;
    int freeBexps;
    TMP_INIT;

    FLINT_ASSERT(k >= 0);

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

    if (A == B)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init(T, ctx);
        nmod_mpoly_pow_si(T, B, k, ctx);
        nmod_mpoly_swap(A, T, ctx);
        nmod_mpoly_clear(T, ctx);
        return;
    }

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

    if (B->length == WORD(1))
    {
        /* powering a monomial */
        nmod_mpoly_fit_length(A, WORD(1), ctx);
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

        fmpz_mpoly_t T;
        slong Tlen;
        fmpz * Bcoeffs_fmpz;

        cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

        T->alloc = k*(B->length - 1) + 1;
        T->coeffs = (fmpz *) flint_calloc(T->alloc, sizeof(fmpz));
        T->exps   = (ulong *) flint_malloc(T->alloc*N*sizeof(ulong));
        T->length = 0;
        T->bits = exp_bits;

        Bcoeffs_fmpz = _fmpz_vec_init(B->length);
        _fmpz_vec_set_nmod_vec(Bcoeffs_fmpz, B->coeffs, B->length, ctx->ffinfo->mod);

        if (!n_is_prime(ctx->ffinfo->mod.n))
        {
            slong Alen, Tlen;

            Tlen = _fmpz_mpoly_pow_fps(&T->coeffs, &T->exps, &T->alloc,
                       Bcoeffs_fmpz, Bexps, B->length, k, exp_bits, N, cmpmask);

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
            ulong ne;
            slong Slen, Ulen;
            nmod_mpoly_t S, R, U;

            nmod_mpoly_init3(S, B->length, exp_bits, ctx);
            nmod_mpoly_init3(U, B->length, exp_bits, ctx);
            nmod_mpoly_init3(R, B->length, exp_bits, ctx);
            nmod_mpoly_one(R, ctx);

            ne = UWORD(1);

            while (k > 0)
            {
                ulong kmodn;
                NMOD_RED(kmodn, (ulong)(k), ctx->ffinfo->mod);

                if (kmodn > 0)
                {
                    /* will accomplish R *= B^(n^e*(k%n)) */

                    if (kmodn <= 2)
                    {
                        nmod_mpoly_fit_length(S, B->length, ctx);

                        /* S = B^(n^e) */
                        for (i = 0; i < B->length; i++)
                        {
                            S->coeffs[i] = B->coeffs[i];
                            mpoly_monomial_mul_ui_mp(S->exps + N*i,
                                                       Bexps + N*i, N, ne);
                        }
                        _nmod_mpoly_set_length(S, B->length, ctx);                        

                        if (kmodn == 2)
                        {
                            /* R *= S */
                            Ulen = _nmod_mpoly_mul_johnson(
                                            &U->coeffs, &U->exps, &U->alloc,
                                             R->coeffs, R->exps, R->length,
                                             S->coeffs, S->exps, S->length,
                                           exp_bits,  N, cmpmask, ctx->ffinfo);
                            _nmod_mpoly_set_length(U, Ulen, ctx);
                            nmod_mpoly_swap(R, U, ctx);
                        }

                    } else
                    {
                        /* S = B^(n^e*(k%n)) */
                        Tlen = _fmpz_mpoly_pow_fps(
                                              &T->coeffs, &T->exps, &T->alloc,
                                              Bcoeffs_fmpz, Bexps, B->length,
                                                  kmodn, exp_bits, N, cmpmask);
                        nmod_mpoly_fit_length(S, Tlen, ctx);
                        Slen = 0;
                        for (i = 0; i < Tlen; i++)
                        {
                            S->coeffs[Slen] = fmpz_fdiv_ui(T->coeffs + i,
                                                           ctx->ffinfo->mod.n);
                            mpoly_monomial_mul_ui_mp(S->exps + N*Slen,
                                                     T->exps + N*i, N, ne);
                            Slen += (S->coeffs[Slen] != UWORD(0));
                        }
                        _nmod_mpoly_set_length(S, Slen, ctx);
                    }

                    if (nmod_mpoly_is_one(R, ctx))
                    {
                        nmod_mpoly_swap(R, S, ctx);
                    } else
                    {
                        Ulen = _nmod_mpoly_mul_johnson(
                                            &U->coeffs, &U->exps, &U->alloc,
                                             R->coeffs, R->exps, R->length,
                                             S->coeffs, S->exps, S->length,
                                           exp_bits,  N, cmpmask, ctx->ffinfo);
                        _nmod_mpoly_set_length(U, Ulen, ctx);
                        nmod_mpoly_swap(R, U, ctx);
                    }
                }

                k = k/ctx->ffinfo->mod.n;
                ne = ne * ctx->ffinfo->mod.n;
            }

            nmod_mpoly_swap(A, R, ctx);

            nmod_mpoly_clear(S, ctx);
            nmod_mpoly_clear(U, ctx);
            nmod_mpoly_clear(R, ctx);
        }

        _fmpz_vec_clear(Bcoeffs_fmpz, B->length);

        for (i = 0; i < T->alloc; i++)
            fmpz_clear(T->coeffs + i);

        flint_free(T->coeffs);
        flint_free(T->exps);

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
