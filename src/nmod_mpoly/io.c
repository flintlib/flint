/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "nmod_mpoly.h"

/* printing *******************************************************************/

static int _nmod_mpoly_fprint_pretty(FILE * file,
                       const mp_limb_t * coeff, const ulong * exp, slong len,
                       const char ** x_in,  slong bits, const mpoly_ctx_t mctx)
{
    slong i, j, N;
    fmpz * exponents;
    int r = 0, first;
    char ** x = (char **) x_in;

    TMP_INIT;

    if (len == 0)
    {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
    }

    N = mpoly_words_per_exp(bits, mctx);

    TMP_START;

    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));
        for (i = 0; i < mctx->nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(((FLINT_BITS+4)/3)*sizeof(char));
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    exponents = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(exponents + i);

    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            r = fputc('+', file);
            r = (r != EOF) ? 1 : EOF;
            if (r <= 0) goto done;
        }

        first = (coeff[i] == 1);
        if (!first)
        {
            r = flint_fprintf(file, "%wu", coeff[i]);
            if (r <= 0) goto done;
        }

        mpoly_get_monomial_ffmpz(exponents, exp + N*i, bits, mctx);

        for (j = 0; j < mctx->nvars; j++)
        {
            int cmp = fmpz_cmp_ui(exponents + j, WORD(1));

            if (cmp < 0)
                continue;

            if (!first)
            {
                r = fputc('*', file);
                r = (r != EOF) ? 1 : EOF;
                if (r <= 0) goto done;
            }

            r = flint_fprintf(file, "%s", x[j]);
            if (r <= 0) goto done;

            if (cmp > 0)
            {
                r = fputc('^', file);
                if (r <= 0) goto done;
                r = fmpz_fprint(file, exponents + j);
                if (r <= 0) goto done;
            }

            first = 0;
        }

        if (first)
        {
            r = flint_fprintf(file, "1");
            if (r <= 0) goto done;
        }
    }

done:
    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;
    return r;
}

int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx) { return _nmod_mpoly_fprint_pretty(file, A->coeffs, A->exps, A->length, x, A->bits, ctx->minfo); }
int nmod_mpoly_print_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx) { return nmod_mpoly_fprint_pretty(stdout, A, x, ctx); }

/* debugging ******************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/

void nmod_mpoly_remainder_strongtest(const nmod_mpoly_t r, const nmod_mpoly_t g,
                                                    const nmod_mpoly_ctx_t ctx)
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
            flint_throw(FLINT_ERROR, "nmod_mpoly_remainder_strongtest FAILED i = %wd\n"
                    "rem %s\n\n"
                    "den %s\n\n",
                    i,
                    nmod_mpoly_get_str_pretty(r, NULL, ctx),
                    nmod_mpoly_get_str_pretty(g, NULL, ctx));
        }
    }

   flint_free(rexp);
   flint_free(gexp);
}

void nmod_mpolyd_print(nmod_mpolyd_t poly)
{

    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (poly->coeffs[i] == 0)
            continue;

        if (!first)
            printf(" + ");

        flint_printf("%wu", poly->coeffs[i]);

        for (j = poly->nvars - 1; j >= 0; j--)
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%wd^%wu", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}

void nmod_mpolyn_print_pretty(const nmod_mpolyn_t A,
                                   const char ** x_in, const nmod_mpoly_ctx_t ctx)
{
    n_poly_struct * coeff = A->coeffs;
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

        printf("(");
        n_poly_print_pretty(coeff + i, "v");
        printf(")");

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
