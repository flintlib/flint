/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpolyn_one(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    n_poly_struct * Acoeff;
    ulong * Aexp;
    slong N;

    fq_nmod_mpolyn_fit_length(A, 1, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);

    n_fq_poly_one(Acoeff + 0, ctx->fqctx);
    mpoly_monomial_zero(Aexp + N*0, N);

    A->length = 1;
}


int fq_nmod_mpolyn_is_canonical(
    const fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
    {
        return 0;
    }

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_poly_is_canonical(A->coeffs + i, ctx->fqctx))
            return 0;

        if (n_poly_is_zero(A->coeffs + i))
            return 0;
    }

    return 1;
}

void fq_nmod_mpolyn_init(fq_nmod_mpolyn_t A, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpolyn_clear(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        n_poly_clear(A->coeffs + i);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_nmod_mpolyn_swap(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B)
{
   fq_nmod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;
}

void fq_nmod_mpolyn_zero(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    A->length = 0;
}

int fq_nmod_mpolyn_is_zero(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length == 0;
}

void fq_nmod_mpolyn_print_pretty(const fq_nmod_mpolyn_t A,
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

    exponents = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(exponents + i);
   
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            printf(" + ");
        }

        flint_printf("(");
        n_fq_poly_print_pretty(A->coeffs + i, "v", ctx->fqctx);
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

void fq_nmod_mpolyn_fit_length(fq_nmod_mpolyn_t A, slong length,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

        A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
        A->coeffs = (n_poly_struct *) flint_realloc(A->coeffs,
                                              new_alloc*sizeof(n_poly_struct));

        for (i = old_alloc; i < new_alloc; i++)
            n_poly_init(A->coeffs + i);

        A->alloc = new_alloc;
    }
}

void fq_nmod_mpolyn_fit_bits(fq_nmod_mpolyn_t A, slong bits, const fq_nmod_mpoly_ctx_t ctx)
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

void fq_nmod_mpolyn_set(fq_nmod_mpolyn_t A, const fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    n_poly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);
/*
    fq_nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;
*/
    Blen = B->length;
    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        n_fq_poly_set(Acoeff + i, Bcoeff + i, ctx->fqctx);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }

    A->length = Blen;
}
