/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

slong _fq_zech_mpoly_sub(fq_zech_struct * coeff1,       ulong * exp1,
                         fq_zech_struct * coeff2, const ulong * exp2, slong len2,
                         fq_zech_struct * coeff3, const ulong * exp3, slong len3,
                   slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx)
{
    slong i = 0, j = 0, k = 0;

    while (i < len2 && j < len3)
    {
        int cmp = mpoly_monomial_cmp(exp2 + i*N, exp3 + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
            fq_zech_set(coeff1 + k, coeff2 + i, fqctx);
            i++;
        } else if (cmp == 0)
        {
            mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
            fq_zech_sub(coeff1 + k, coeff2 + i, coeff3 + j, fqctx);
            k -= fq_zech_is_zero(coeff1 + k, fqctx);
            i++;
            j++;
        } else
        {
            mpoly_monomial_set(exp1 + k*N, exp3 + j*N, N);
            fq_zech_neg(coeff1 + k, coeff3 + j, fqctx);
            j++;
        }
        k++;
    }

    while (i < len2)
    {
        mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
        fq_zech_set(coeff1 + k, coeff2 + i, fqctx);
        i++;
        k++;
    }

    while (j < len3)
    {
        mpoly_monomial_set(exp1 + k*N, exp3 + j*N, N);
        fq_zech_neg(coeff1 + k, coeff3 + j, fqctx);
        j++;
        k++;
    }

    return k;
}

void fq_zech_mpoly_sub(fq_zech_mpoly_t poly1, const fq_zech_mpoly_t poly2,
                          const fq_zech_mpoly_t poly3, const fq_zech_mpoly_ctx_t ctx)
{
    slong len1 = 0, max_bits, N;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    ulong * cmpmask;
    int free2 = 0, free3 = 0;
    TMP_INIT;

    max_bits = FLINT_MAX(poly2->bits, poly3->bits);
    N = mpoly_words_per_exp(max_bits, ctx->minfo);

    if (poly2->length == 0)
    {
        fq_zech_mpoly_neg(poly1, poly3, ctx);
        return;
    }
    else if (poly3->length == 0)
    {
        fq_zech_mpoly_set(poly1, poly2, ctx);
        return;
    }

    TMP_START;
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, max_bits, ctx->minfo);

    if (max_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, max_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
    }

    if (max_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, max_bits, poly3->exps, poly3->bits,
                                                    poly3->length, ctx->minfo);
    }

    if (poly1 == poly2 || poly1 == poly3)
    {
        fq_zech_mpoly_t temp;

        fq_zech_mpoly_init2(temp, poly2->length + poly3->length, ctx);
        fq_zech_mpoly_fit_bits(temp, max_bits, ctx);
        temp->bits = max_bits;

        len1 = _fq_zech_mpoly_sub(temp->coeffs, temp->exps, 
                    poly2->coeffs, exp2, poly2->length,
                    poly3->coeffs, exp3, poly3->length,
                                    N, cmpmask, ctx->fqctx);

        fq_zech_mpoly_swap(temp, poly1, ctx);
        fq_zech_mpoly_clear(temp, ctx);

    }
    else
    {
        fq_zech_mpoly_fit_length(poly1, poly2->length + poly3->length, ctx);
        fq_zech_mpoly_fit_bits(poly1, max_bits, ctx);
        poly1->bits = max_bits;

        len1 = _fq_zech_mpoly_sub(poly1->coeffs, poly1->exps, 
                       poly2->coeffs, exp2, poly2->length,
                       poly3->coeffs, exp3, poly3->length,
                                    N, cmpmask, ctx->fqctx);
    }
      
    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    _fq_zech_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
