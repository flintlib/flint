/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fq_nmod_mpoly.h"

/* printing *******************************************************************/

int fq_nmod_mpoly_fprint_pretty(
    FILE * file,
    const fq_nmod_mpoly_t A,
    const char ** x_in,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong len = A->length;
    ulong * exp = A->exps;
    slong bits = A->bits;
    slong i, j, N;
    fmpz * exponents;
    int r = 0;
    char ** x = (char **) x_in;
    TMP_INIT;

    if (len == 0)
    {
        r = (EOF != fputc('0', file));
        return r;
    }

#define CHECK_r if (r <= 0) goto done;

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
            r = flint_fprintf(file, " + ");
            CHECK_r
        }

        r = flint_fprintf(file, "(");
        CHECK_r
        r = n_fq_fprint_pretty(file, A->coeffs + d*i, ctx->fqctx);
        CHECK_r
        r = flint_fprintf(file, ")");
        CHECK_r

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, ctx->minfo);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, 1);

            if (cmp > 0)
            {
                r = flint_fprintf(file, "*%s^", x[j]);
                CHECK_r
                r = fmpz_fprint(file, exponents + j);
                CHECK_r
            }
            else if (cmp == 0)
            {
                r = flint_fprintf(file, "*%s", x[j]);
                CHECK_r
            }
        }
    }

done:

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    return r;
}

int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx) { return fq_nmod_mpoly_fprint_pretty(stdout, A, x, ctx); }

/* debugging ******************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/

void fq_nmod_mpoly_remainder_strongtest(const fq_nmod_mpoly_t r, const fq_nmod_mpoly_t g, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N, bits;
    ulong mask = 0;
    ulong * rexp, * gexp;

    bits = FLINT_MAX(r->bits, g->bits);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    if (g->length == 0 )
        flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

    if (r->length == 0 )
        return;

    rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
    gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
    mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
    mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

    if (bits <= FLINT_BITS)
        mask = mpoly_overflow_mask_sp(bits);
    else
        mask = 0;

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_throw(FLINT_ERROR, "fq_nmod_mpoly_remainder_strongtest FAILED i = %wd\n"
                    "rem %s\n\n"
                    "den %s\n\n",
                    i,
                    fq_nmod_mpoly_get_str_pretty(r, NULL, ctx),
                    fq_nmod_mpoly_get_str_pretty(g, NULL, ctx));
        }
    }

    flint_free(rexp);
    flint_free(gexp);
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
