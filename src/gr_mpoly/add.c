/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

static slong _gr_mpoly_add(
    slong * Alen,
    gr_ptr Acoeffs, ulong * Aexps,
    gr_srcptr Bcoeffs, const ulong * Bexps, slong Blen,
    gr_srcptr Ccoeffs, const ulong * Cexps, slong Clen,
    slong N,
    const ulong * cmpmask,
    gr_ctx_t fctx)
{
    gr_method_binary_op add = GR_BINARY_OP(fctx, ADD);
    gr_method_unary_op set = GR_UNARY_OP(fctx, SET);
    slong sz = fctx->sizeof_elem;
    slong i = 0, j = 0, k = 0;

    int status = GR_SUCCESS;

    while (i < Blen && j < Clen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + i*N, Cexps + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            status |= set(GR_ENTRY(Acoeffs, k, sz), GR_ENTRY(Bcoeffs, i, sz), fctx);
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            status |= add(GR_ENTRY(Acoeffs, k, sz), GR_ENTRY(Bcoeffs, i, sz), GR_ENTRY(Ccoeffs, j, sz), fctx);
            k += (gr_is_zero(GR_ENTRY(Acoeffs, k, sz), fctx) != T_TRUE);
            i++;
            j++;
        }
        else
        {
            mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
            status |= set(GR_ENTRY(Acoeffs, k, sz), GR_ENTRY(Ccoeffs, j, sz), fctx);
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        status |= set(GR_ENTRY(Acoeffs, k, sz), GR_ENTRY(Bcoeffs, i, sz), fctx);
        i++;
        k++;
    }

    while (j < Clen)
    {
        mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
        status |= set(GR_ENTRY(Acoeffs, k, sz), GR_ENTRY(Ccoeffs, j, sz), fctx);
        j++;
        k++;
    }

    *Alen = k;

    return status;
}

int gr_mpoly_add(
    gr_mpoly_t A,
    const gr_mpoly_t B,
    const gr_mpoly_t C,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong N;
    flint_bitcnt_t Abits;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    ulong * cmpmask;
    int freeBexps = 0, freeCexps = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    if (B->length == 0)
        return gr_mpoly_set(A, C, mctx, cctx);

    if (C->length == 0)
        return gr_mpoly_set(A, B, mctx, cctx);

    TMP_START;
    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, mctx);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, mctx);

    if (Abits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits,
                                                    B->length, mctx);
    }

    if (Abits != C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits,
                                                    C->length, mctx);
    }

    if (A == B || A == C)
    {
        gr_mpoly_t T;
        gr_mpoly_init3(T, B->length + C->length, Abits, mctx, cctx);
        status = _gr_mpoly_add(&T->length, T->coeffs, T->exps,
                                        B->coeffs, Bexps, B->length,
                                        C->coeffs, Cexps, C->length,
                                                      N, cmpmask, cctx);
        gr_mpoly_swap(A, T, mctx, cctx);
        gr_mpoly_clear(T, mctx, cctx);
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, mctx, cctx);
        status = _gr_mpoly_add(&A->length, A->coeffs, A->exps,
                                        B->coeffs, Bexps, B->length,
                                        C->coeffs, Cexps, C->length,
                                                      N, cmpmask, cctx);
    }

    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;

    return status;
}
