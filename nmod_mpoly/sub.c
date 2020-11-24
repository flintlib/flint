/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_sub1(mp_limb_t * coeff1,       ulong * exp1,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
                                          ulong maskhi, nmod_t fctx)
{
    slong i = 0, j = 0, k = 0;

    while (i < len2 && j < len3)
    {
        if ((exp2[i]^maskhi) > (exp3[j]^maskhi))
        {
            exp1[k] = exp2[i];
            coeff1[k] = coeff2[i];
            i++;
        } else if ((exp2[i]^maskhi) == (exp3[j]^maskhi))
        {
            exp1[k] = exp2[i];
            coeff1[k] = nmod_sub(coeff2[i], coeff3[j], fctx);
            k -= (coeff1[k] == 0);
            i++;
            j++;
        } else
        {
            exp1[k] = exp3[j];
            coeff1[k] = nmod_neg(coeff3[j], fctx);
            j++;         
        }
        k++;
    }

    while (i < len2)
    {
        exp1[k] = exp2[i];
        coeff1[k] = coeff2[i];
        i++;
        k++;
    }

    while (j < len3)
    {
        exp1[k] = exp3[j];
        coeff1[k] = nmod_neg(coeff3[j], fctx);
        j++;
        k++;
    }

    return k;
}

slong _nmod_mpoly_sub(ulong * coeff1,       ulong * exp1,
                const ulong * coeff2, const ulong * exp2, slong len2,
                const ulong * coeff3, const ulong * exp3, slong len3,
                   slong N, const ulong * cmpmask, nmod_t fctx)
{
    slong i = 0, j = 0, k = 0;

    if (N == 1)
        return _nmod_mpoly_sub1(coeff1, exp1, coeff2, exp2, len2,
                                         coeff3, exp3, len3, cmpmask[0], fctx);

    while (i < len2 && j < len3)
    {
        int cmp = mpoly_monomial_cmp(exp2 + i*N, exp3 + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
            coeff1[k] = coeff2[i];
            i++;
        } else if (cmp == 0)
        {
            mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
            coeff1[k] = nmod_sub(coeff2[i], coeff3[j], fctx);
            k -= (coeff1[k] == 0);
            i++;
            j++;
        } else
        {
            mpoly_monomial_set(exp1 + k*N, exp3 + j*N, N);
            coeff1[k] = nmod_neg(coeff3[j], fctx);
            j++;
        }
        k++;
    }

    while (i < len2)
    {
        mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
        coeff1[k] = coeff2[i];
        i++;
        k++;
    }

    while (j < len3)
    {
        mpoly_monomial_set(exp1 + k*N, exp3 + j*N, N);
        coeff1[k] = nmod_neg(coeff3[j], fctx);
        j++;
        k++;
    }

    return k;
}

void nmod_mpoly_sub(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                          const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
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
        nmod_mpoly_neg(poly1, poly3, ctx);
        return;
    } else if (poly3->length == 0)
    {
        nmod_mpoly_set(poly1, poly2, ctx);
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
        nmod_mpoly_t temp;
        nmod_mpoly_init3(temp, poly2->length + poly3->length, max_bits, ctx);
        len1 = _nmod_mpoly_sub(temp->coeffs, temp->exps, 
                                    poly2->coeffs, exp2, poly2->length,
                                    poly3->coeffs, exp3, poly3->length,
                                                      N, cmpmask, ctx->mod);
        nmod_mpoly_swap(temp, poly1, ctx);
        nmod_mpoly_clear(temp, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(poly1, poly2->length + poly3->length, max_bits, ctx);
        len1 = _nmod_mpoly_sub(poly1->coeffs, poly1->exps, 
                                    poly2->coeffs, exp2, poly2->length,
                                    poly3->coeffs, exp3, poly3->length,
                                                      N, cmpmask, ctx->mod);
    }
      
    _nmod_mpoly_set_length(poly1, len1, ctx);

    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    TMP_END;
}
