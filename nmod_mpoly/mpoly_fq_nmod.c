/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t poly, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        return 0;

    if (mpoly_monomials_overflow_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        return 0;

    for (i = 0; i < poly->length; i++)
    {
        if (fq_nmod_is_zero(poly->coeffs + i, ctx->fqctx))
            return 0;
    }

    return 1;
}


void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                                        mp_limb_t p, slong deg)
{
    fmpz_t P;

    mpoly_ctx_init(ctx->minfo, nvars, ORD_LEX);

    fmpz_init_set_ui(P, p);
    fq_nmod_ctx_init(ctx->fqctx, P, deg, "#");
    fmpz_clear(P);
}

void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_ctx_clear(ctx->minfo);
    fq_nmod_ctx_clear(ctx->fqctx);
}

void fq_nmod_mpoly_ctx_change_modulus(fq_nmod_mpoly_ctx_t ctx, slong deg)
{
    fmpz_t P;
    fmpz_init_set_ui(P, ctx->fqctx->mod.n);
    fq_nmod_ctx_clear(ctx->fqctx);
    fq_nmod_ctx_init(ctx->fqctx, P, deg, "#");
    fmpz_clear(P);
}

void fq_nmod_mpoly_init(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong bits = mpoly_fix_bits(WORD(8), ctx->minfo);

    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpoly_init2(fq_nmod_mpoly_t poly,
                                       slong alloc, const fq_nmod_mpoly_ctx_t ctx)
{
    /* default to at least 8 bits per exponent */
    slong bits = mpoly_fix_bits(WORD(8), ctx->minfo);
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (alloc != 0)
    {
        slong i;
        poly->coeffs = (fq_nmod_struct *) flint_malloc(alloc*sizeof(fq_nmod_struct));
        poly->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
        for (i = 0; i < alloc; i++)
            fq_nmod_init(poly->coeffs + i, ctx->fqctx);
    } else
    {
        poly->coeffs = NULL;
        poly->exps = NULL;
    }
    poly->alloc = alloc;
    poly->length = 0;
    poly->bits = bits;
}

void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_clear(A->coeffs + i, ctx->fqctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_nmod_clear(A->coeffs + i, ctx->fqctx);
        fq_nmod_init(A->coeffs + i, ctx->fqctx);
    }
    A->length = 0;
}

int fq_nmod_mpoly_is_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length == 0;
}

int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t poly, const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (poly->length != 1)
        return 0;

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(poly->exps + N*0, N))
        return 0;

    return fq_nmod_is_one(poly->coeffs + 0, ctx->fqctx);
}


void fq_nmod_mpoly_swap(fq_nmod_mpoly_t poly1, fq_nmod_mpoly_t poly2)
{
   fq_nmod_mpoly_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

void _fq_nmod_mpoly_set_length(fq_nmod_mpoly_t poly, slong newlen, 
                                                   const fq_nmod_mpoly_ctx_t ctx)
{
    poly->length = newlen;
}

void fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A,
                             const char ** x_in, const fq_nmod_mpoly_ctx_t ctx)
{
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        flint_printf("0");
        return;
    }

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(ctx->minfo->nvars*sizeof(char *));
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i+1);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            printf(" + ");
        }

        flint_printf("(");
        fq_nmod_print_pretty(A->coeffs + i, ctx->fqctx);
        flint_printf(")");

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));

            if (cmp > 0)
            {
                printf("*%s^", x[j]);
                fmpz_print(exponents + j);
            } else if (cmp == 0)
            {
                printf("*%s", x[j]);
            }
        }
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
}

void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            A->coeffs = (fq_nmod_struct *) flint_malloc(new_alloc*sizeof(fq_nmod_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            A->coeffs = (fq_nmod_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fq_nmod_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fq_nmod_init(A->coeffs + i, ctx->fqctx);
        }
        A->alloc = new_alloc;
    }
}

void _fq_nmod_mpoly_fit_length(fq_nmod_struct ** coeff,
                              ulong ** exps, slong * alloc, slong len, slong N, const fq_nmod_ctx_t fqctx)
{
    if (len > *alloc)
    {
        slong i;
        len = FLINT_MAX(len, 2*(*alloc));
        (* coeff) = (fq_nmod_struct *) flint_realloc(* coeff, len*sizeof(fq_nmod_struct));
        (* exps) = (ulong *) flint_realloc(*exps, len*N*sizeof(ulong)); 
        for (i = *alloc; i < len; i++)
            fq_nmod_init((* coeff) + i, fqctx);
        (* alloc) = len;
    }
}



void fq_nmod_mpoly_fit_bits(fq_nmod_mpoly_t A, slong bits, const fq_nmod_mpoly_ctx_t ctx)
{
   slong N;
   ulong * t;

   if (A->bits < bits)
   {
      if (A->alloc != 0)
      {
         N = mpoly_words_per_exp(bits, ctx->minfo);
         t = flint_malloc(N*A->alloc*sizeof(ulong));
         mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, ctx->minfo);
         flint_free(A->exps);
         A->exps = t;
      }

      A->bits = bits;
   }
}

void fq_nmod_mpoly_set_length(fq_nmod_mpoly_t A, slong newlen, const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;
        for (i = newlen; i < A->length; i++)
        {
            fq_nmod_clear(A->coeffs + i, ctx->fqctx);
            fq_nmod_init(A->coeffs + i, ctx->fqctx);
        }
    }
    A->length = newlen;
}





void _fq_nmod_mpoly_set(fq_nmod_struct * coeff1, ulong * exps1,
                const fq_nmod_struct * coeff2, const ulong * exps2, slong len2,
                                               slong N, const fq_nmod_ctx_t fqctx)
{
    slong i;

    if (coeff1 != coeff2)
    {
        for (i = 0; i < len2; i++)
            fq_nmod_set(coeff1 + i, coeff2 + i, fqctx);
    }

    if (exps1 != exps2)
    {
        for (i = 0; i < len2; i++)
            mpoly_monomial_set(exps1 + N*i, exps2 + N*i, N);
    }
}

void fq_nmod_mpoly_set(fq_nmod_mpoly_t poly1, const fq_nmod_mpoly_t poly2,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    fq_nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);

    _fq_nmod_mpoly_set(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length, N, ctx->fqctx);

    _fq_nmod_mpoly_set_length(poly1, poly2->length, ctx);
    poly1->bits = poly2->bits;
}



