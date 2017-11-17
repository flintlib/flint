/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_add1(ulong * coeff1,       ulong * exp1,
                 const ulong * coeff2, const ulong * exp2, slong len2,
                 const ulong * coeff3, const ulong * exp3, slong len3,
                                          ulong maskhi, const nmodf_ctx_t fctx)
{
    slong i = 0, j = 0, k = 0;
    slong D = fctx->deg;

    while (i < len2 && j < len3)
    {
        if ((exp2[i]^maskhi) > (exp3[j]^maskhi))
        {
            exp1[k] = exp2[i];
            nmodf_set(coeff1 + k*D, coeff2 + i*D, fctx);
            i++;
        } else if ((exp2[i]^maskhi) == (exp3[j]^maskhi))
        {
            exp1[k] = exp2[i];
            nmodf_add(coeff1 + k*D, coeff2 + i*D, coeff3 + j*D, fctx);
            if (nmodf_is_zero(coeff1 + k*D, fctx))
                k--;
            i++;
            j++;
        } else
        {
            nmodf_set(coeff1 + k*D, coeff3 + j*D, fctx);
            exp1[k] = exp3[j];
            j++;         
        }
        k++;
    }

    while (i < len2)
    {
        exp1[k] = exp2[i];
        nmodf_set(coeff1 + k*D, coeff2 + i*D, fctx);
        i++;
        k++;
    }

    while (j < len3)
    {
        exp1[k] = exp3[j];
        nmodf_set(coeff1 + k*D, coeff3 + j*D, fctx);
        j++;
        k++;
    }

    return k;
}

slong _nmod_mpoly_add(ulong * coeff1,       ulong * exp1,
                const ulong * coeff2, const ulong * exp2, slong len2,
                const ulong * coeff3, const ulong * exp3, slong len3,
                   slong N, ulong maskhi, ulong masklo, const nmodf_ctx_t fctx)
{
    slong i = 0, j = 0, k = 0;
    slong D = fctx->deg;

    if (N == 1)
        return _nmod_mpoly_add1(coeff1, exp1, coeff2, exp2, len2,
                                            coeff3, exp3, len3, maskhi, fctx);

    while (i < len2 && j < len3)
    {
        int cmp = mpoly_monomial_cmp(exp2 + i*N, exp3 + j*N, N, maskhi, masklo);

        if (cmp > 0)
        {
            mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
            nmodf_set(coeff1 + k*D, coeff2 + i*D, fctx);
            i++;
        } else if (cmp == 0)
        {
            mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
            nmodf_add(coeff1 + k*D, coeff2 + i*D, coeff3 + j*D, fctx);
            if (nmodf_is_zero(coeff1 + k*D, fctx))
                k--;
            i++;
            j++;
        } else
        {
            nmodf_set(coeff1 + k*D, coeff3 + j*D, fctx);
            mpoly_monomial_set(exp1 + k*N, exp3 + j*N, N);
            j++;         
        }
        k++;
    }

    while (i < len2)
    {
        mpoly_monomial_set(exp1 + k*N, exp2 + i*N, N);
        nmodf_set(coeff1 + k*D, coeff2 + i*D, fctx);
        i++;
        k++;
    }

    while (j < len3)
    {
        mpoly_monomial_set(exp1 + k*N, exp3 + j*N, N);
        nmodf_set(coeff1 + k*D, coeff3 + j*D, fctx);
        j++;
        k++;
    }

    return k;
}

void nmod_mpoly_add(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                          const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
{
    slong len1 = 0, max_bits, N;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    ulong maskhi, masklo;
    int free2 = 0, free3 = 0;

    max_bits = FLINT_MAX(poly2->bits, poly3->bits);
    masks_from_bits_ord(maskhi, masklo, max_bits, ctx->ord);
    N = words_per_exp(ctx->n, max_bits);

    if (poly2->length == 0)
    {
        nmod_mpoly_set(poly1, poly3, ctx);
        return;
    } else if (poly3->length == 0)
    {
        nmod_mpoly_set(poly1, poly2, ctx);
        return;
    }

    if (max_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_unpack_monomials(exp2, max_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
    }

    if (max_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_unpack_monomials(exp3, max_bits, poly3->exps, poly3->bits,
                                                        poly3->length, ctx->n);
    }

    if (poly1 == poly2 || poly1 == poly3)
    {
        nmod_mpoly_t temp;

        nmod_mpoly_init2(temp, poly2->length + poly3->length, ctx);
        nmod_mpoly_fit_bits(temp, max_bits, ctx);
        temp->bits = max_bits;

        len1 = _nmod_mpoly_add(temp->coeffs, temp->exps, 
                    poly2->coeffs, exp2, poly2->length,
                    poly3->coeffs, exp3, poly3->length,
                                    N, maskhi, masklo, ctx->ffinfo);

        nmod_mpoly_swap(temp, poly1, ctx);

        nmod_mpoly_clear(temp, ctx);
    } else
    {
        nmod_mpoly_fit_length(poly1, poly2->length + poly3->length, ctx);
        nmod_mpoly_fit_bits(poly1, max_bits, ctx);
        poly1->bits = max_bits;

        len1 = _nmod_mpoly_add(poly1->coeffs, poly1->exps, 
                       poly2->coeffs, exp2, poly2->length,
                       poly3->coeffs, exp3, poly3->length,
                                    N, maskhi, masklo, ctx->ffinfo);
    }
      
    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    _nmod_mpoly_set_length(poly1, len1, ctx);
}
