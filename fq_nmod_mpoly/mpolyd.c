/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpolyd_ctx_init(fq_nmod_mpolyd_ctx_t dctx, slong nvars,
                                                        mp_limb_t p, slong deg)
{
    slong i;
    fmpz_t P;
    fmpz_init_set_ui(P, p);

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    fq_nmod_ctx_init(dctx->fqctx, P, deg, "#");
    fmpz_clear(P);
}

void fq_nmod_mpolyd_ctx_init_modulus(fq_nmod_mpolyd_ctx_t dctx, slong nvars,
                                                    const fq_nmod_ctx_t fqctx)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    fq_nmod_ctx_init_modulus(dctx->fqctx, fqctx->modulus, "#");
}

void fq_nmod_mpolyd_ctx_init2(fq_nmod_mpolyd_ctx_t dctx, slong nvars,
                                                     const fq_nmod_ctx_t fqctx)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    fq_nmod_ctx_init_modulus(dctx->fqctx, fqctx->modulus, fqctx->var);
}

void fq_nmod_mpolyd_ctx_clear(fq_nmod_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
    fq_nmod_ctx_clear(dctx->fqctx);
}

void fq_nmod_mpolyd_init(fq_nmod_mpolyd_t poly, slong nvars,
                                                     const fq_nmod_ctx_t fqctx)
{
    slong i;

    poly->nvars = nvars;
    poly->degb_alloc = nvars;
    poly->deg_bounds = (slong *) flint_malloc(poly->degb_alloc*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeff_alloc = WORD(16);
    poly->coeffs = (fq_nmod_struct *) flint_malloc(poly->coeff_alloc
                                                      *sizeof(fq_nmod_struct));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fq_nmod_init(poly->coeffs + i, fqctx);
    }
}

void fq_nmod_mpolyd_fit_length(fq_nmod_mpolyd_t poly, slong len,
                                                     const fq_nmod_ctx_t fqctx)
{
    slong i;
    if (poly->coeff_alloc < len) {
        slong new_alloc = len;
        poly->coeffs = (fq_nmod_struct *) flint_realloc(poly->coeffs,
                                             new_alloc*sizeof(fq_nmod_struct));
        for (i = poly->coeff_alloc; i < new_alloc; i++)
        {
            fq_nmod_init(poly->coeffs + i, fqctx);
        }        
        poly->coeff_alloc = new_alloc;
    }
}

void fq_nmod_mpolyd_set_nvars(fq_nmod_mpolyd_t poly, slong nvars)
{
    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars
                                                               *sizeof(slong));
        poly->degb_alloc = nvars;
    }
}

void fq_nmod_mpolyd_zero(fq_nmod_mpolyd_t poly, const fq_nmod_ctx_t fqctx)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    fq_nmod_zero(poly->coeffs + 0, fqctx);
}

void fq_nmod_mpolyd_clear(fq_nmod_mpolyd_t poly, const fq_nmod_ctx_t fqctx)
{
    slong i;
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fq_nmod_clear(poly->coeffs + i, fqctx);
    }
    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}

void fq_nmod_mpolyd_print(fq_nmod_mpolyd_t poly, const fq_nmod_ctx_t fqctx)
{
    int first = 0;
    slong i, j;
    slong degb_prod;

    flint_printf("[ ");
    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++)
    {
        flint_printf("%wd ", poly->deg_bounds[j]);
        degb_prod *= poly->deg_bounds[j];
    }
    flint_printf("]: ");

    first = 1;
    for (i = 0; i < degb_prod; i++)
    {
        ulong k = i;

        if (fq_nmod_is_zero(poly->coeffs + i, fqctx))
            continue;

        if (!first)
            printf(" + ");

        flint_printf("(");
        fq_nmod_print_pretty(poly->coeffs + i, fqctx);
        flint_printf(")");

        for (j = poly->nvars - 1; j >= 0; j--) 
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%wd^%wd", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
    {
        flint_printf("0");
    }
}


/*
    m is the number of variables in A
*/
void fq_nmod_mpoly_to_fq_nmod_mpolyd_perm_deflate(fq_nmod_mpolyd_t A, slong m,
              const fq_nmod_mpoly_t B, const slong * perm, const ulong * shift,
     const ulong * stride, const ulong * degree, const fq_nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong degb_prod;
    slong i, k, l, N;
    ulong * Bexp;
    TMP_INIT;

    FLINT_ASSERT(m <= n);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(B->length > 0);

    fq_nmod_mpolyd_set_nvars(A, m);

    TMP_START;
    Bexp = (ulong *) TMP_ALLOC(n*sizeof(slong));

    degb_prod = WORD(1);
    for (k = 0; k < m; k++)
    {
        l = perm[k];
        FLINT_ASSERT(stride[l] != UWORD(0));
        FLINT_ASSERT((degree[l] - shift[l]) % stride[l] == UWORD(0));
        A->deg_bounds[k] = (degree[l] - shift[l])/stride[l] + 1;
        degb_prod *= A->deg_bounds[k];
        /* we should not be converting something whose dense size overflows */
        FLINT_ASSERT(degb_prod > 0);
        FLINT_ASSERT(degb_prod >= A->deg_bounds[k]);
    }

    fq_nmod_mpolyd_fit_length(A, degb_prod, ctx->fqctx);
    for (i = 0; i < degb_prod; i++)
    {
        fq_nmod_zero(A->coeffs + i, ctx->fqctx);
    }

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        slong off = 0;
        mpoly_get_monomial_ui(Bexp, B->exps + N*i, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexp[l] - shift[l]) % stride[l]) == UWORD(0));
            FLINT_ASSERT((Bexp[l] - shift[l])/stride[l] < A->deg_bounds[k]);
            off = (Bexp[l] - shift[l])/stride[l] + A->deg_bounds[k]*off;
        }
        fq_nmod_set(A->coeffs + off, B->coeffs + i, ctx->fqctx);
    }

    TMP_END;
}


void fq_nmod_mpoly_from_fq_nmod_mpolyd_perm_inflate(fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx, const fq_nmod_mpolyd_t B,
                 const slong * perm, const ulong * shift, const ulong * stride)
{
    slong off;
    slong n = ctx->minfo->nvars;
    slong m = B->nvars;
    slong Alen;
    slong i, j, l, k, N;
    slong perm_nontrivial;
    ulong topmask;
    ulong * exps, * pcurexp, * pexps;
    TMP_INIT;

    FLINT_ASSERT(m <= n);
    FLINT_ASSERT(Abits <= FLINT_BITS);

    perm_nontrivial = n - m;

    /* we are going to push back terms manually */
    Alen = 0;
    fq_nmod_mpoly_zero(A, ctx);
    fq_nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    N = mpoly_words_per_exp(Abits, ctx->minfo);

    TMP_START;

    /* find exponent vector for all variables in B */
    pexps = (ulong *) TMP_ALLOC(N*m*sizeof(ulong));
    for (k = 0; k < m; k++)
    {
        l = perm[k];
        perm_nontrivial |= l - k;
        mpoly_gen_monomial_sp(pexps + k*N, l, Abits, ctx->minfo);
        mpoly_monomial_mul_ui(pexps + k*N, pexps + k*N, N, stride[l]);
    }

    /* get most significant exponent in pcurexp and its vector in exps */
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    off = WORD(1);
    for (j = 0; j < m; j++)
    {
        off *= B->deg_bounds[j];
    }
    FLINT_ASSERT(off <= B->coeff_alloc);
    off--;
    mpoly_set_monomial_ui(pcurexp, shift, Abits, ctx->minfo);
    i = off;
    for (k = m - 1; k >= 0; k--) 
    {
        exps[k] = i % B->deg_bounds[k];
        i = i / B->deg_bounds[k];
        mpoly_monomial_madd(pcurexp, pcurexp, exps[k], pexps + N*k, N);
    }

    /* scan down through the exponents */
    topmask = 0;

    for (; off >= 0; off--)
    {
        if (!fq_nmod_is_zero(B->coeffs + off, ctx->fqctx))
        {
            _fq_nmod_mpoly_fit_length(&A->coeffs, &A->exps, &A->alloc,
                                                      Alen + 1, N, ctx->fqctx);
            fq_nmod_set(A->coeffs + Alen, B->coeffs + off, ctx->fqctx);
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= pcurexp[N - 1];
            Alen++;
        }

        k = m - 1;
        do {
            --exps[k];
            if ((slong)(exps[k]) < WORD(0))
            {
                FLINT_ASSERT(off == 0 || k > 0);
                FLINT_ASSERT(exps[k] == -UWORD(1));
                exps[k] = B->deg_bounds[k] - 1;
                mpoly_monomial_madd(pcurexp, pcurexp, exps[k], pexps + N*k, N);
            }
            else
            {
                mpoly_monomial_sub(pcurexp, pcurexp, pexps + N*k, N);
                break;
            }
        } while (--k >= 0);
    }
    _fq_nmod_mpoly_set_length(A, Alen, ctx);

    /* sort the terms if needed */
    if (ctx->minfo->ord != ORD_LEX || perm_nontrivial != WORD(0))
    {
        slong msb;
        mpoly_get_cmpmask(pcurexp, N, Abits, ctx->minfo);
        if (topmask != WORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        }
        else
        {
            msb = -WORD(1);
        }
        if (N == 1)
        {
            if (msb >= WORD(0))
            {
                _fq_nmod_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, pcurexp[0], topmask);
            }
        }
        else
        {
            _fq_nmod_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, pcurexp);
        }
    }

    TMP_END;
}


/*
    convert B to A - sets degree bounds in A
*/
void fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(
                          fq_nmod_mpolyd_t A, const fq_nmod_mpolyd_ctx_t dctx,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    fq_nmod_mpolyd_set_nvars(A, ctx->minfo->nvars);

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (B->length == 0)
    {
        fq_nmod_mpolyd_zero(A, dctx->fqctx);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    fq_nmod_mpoly_degrees_si(exps, B, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        A->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= A->deg_bounds[i];
    }

    fq_nmod_mpolyd_fit_length(A, degb_prod, dctx->fqctx);
    for (i = 0; i < degb_prod; i++)
    {
        fq_nmod_zero(A->coeffs + i, dctx->fqctx);
    }

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, B->exps + N*i, B->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + A->deg_bounds[j]*off;
        }

        fq_nmod_set(A->coeffs + off, B->coeffs + i, ctx->fqctx);
    }

    TMP_END;
}

/*
    Convert B to A
*/
void fq_nmod_mpoly_convert_from_fq_nmod_mpolyd(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx,
                           const fq_nmod_mpolyd_t B, const fq_nmod_mpolyd_ctx_t dctx)
{
    slong off, j, k, N;
    slong bits, nvars = ctx->minfo->nvars;
    slong Alen;
    slong * perm = dctx->perm;
    slong perm_nontrivial = 0;
    ulong topmask;
    ulong * exps, * pcurexp, * pexps;
    TMP_INIT;

    FLINT_ASSERT(nvars == B->nvars);

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    /* find bits needed for the result */
    off = 1;
    for (j = 0; j < nvars; j++)
    {
        off *= B->deg_bounds[j];
        exps[perm[j]] = B->deg_bounds[j] - 1;
        perm_nontrivial |= j ^ perm[j];
    }

    FLINT_ASSERT(off <= B->coeff_alloc);

    bits = mpoly_exp_bits_required_ui(exps, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* we are going to push back terms manually */
    Alen = 0;
    fq_nmod_mpoly_zero(A, ctx);
    fq_nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    /* find exponent vector for all variables */
    pexps = (ulong *) TMP_ALLOC(N*nvars*sizeof(ulong));
    for (k = 0; k < nvars; k++)
    {
        for (j = 0; j < nvars; j++)
            exps[perm[j]] = (j == k);
        mpoly_set_monomial_ui(pexps + k*N, exps, bits, ctx->minfo);
    }

    /* get most significant exponent in exps and its vector in ptempexp */
    off--;
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(pcurexp, N);
    k = off;
    for (j = nvars - 1; j >= 0; j--) 
    {
        exps[j] = k % B->deg_bounds[j];
        k = k / B->deg_bounds[j];
        mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
    }

    /* scan down through the exponents */
    topmask = 0;
    for (; off >= 0; off--)
    {
        if (!fq_nmod_is_zero(B->coeffs + off, ctx->fqctx))
        {
            _fq_nmod_mpoly_fit_length(&A->coeffs, &A->exps, &A->alloc, Alen + 1, N, ctx->fqctx);
            fq_nmod_set(A->coeffs + Alen, B->coeffs + off, dctx->fqctx);
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= (A->exps + N*Alen)[N - 1];
            Alen++;
        }

        j = nvars - 1;
        do {
            --exps[j];
            if ((slong)(exps[j]) < WORD(0))
            {
                FLINT_ASSERT(off == 0 || j > 0);
                FLINT_ASSERT(exps[j] == -UWORD(1));
                exps[j] = B->deg_bounds[j] - 1;
                mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
            }
            else
            {
                mpoly_monomial_sub_mp(pcurexp, pcurexp, pexps + N*j, N);
                break;
            }
        } while (--j >= 0);
    }
    _fq_nmod_mpoly_set_length(A, Alen, ctx);

    /* sort the exponents if needed */
    if (ctx->minfo->ord != ORD_LEX || perm_nontrivial != WORD(0))
    {
        slong msb;
        mpoly_get_cmpmask(pcurexp, N, bits, ctx->minfo);
        if (topmask != WORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        }
        else
        {
            msb = -WORD(1);
        }
        if (N == 1)
        {
            if (msb >= WORD(0))
            {
                _fq_nmod_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, pcurexp[0], topmask);
            }
        }
        else
        {
            _fq_nmod_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, pcurexp);
        }
    }

    TMP_END;
}
