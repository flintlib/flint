/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "mpoly.h"
#include "fmpz_mod_mpoly_factor.h"

/* printing *******************************************************************/

void fmpz_mod_polyu2n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 1, 2),
            var1, extract_exp(A->exps[i], 0, 2));
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_polyu1n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * varlast,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu", var0, A->exps[i]);
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_polyu3n_print_pretty(
    const fmpz_mod_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_mpoly_factor_print_pretty(
    const fmpz_mod_mpoly_factor_t f,
    const char ** vars,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    fmpz_print(f->constant);
    for (i = 0; i < f->num; i++)
    {
        flint_printf("\n*(", i);
        fmpz_mod_mpoly_print_pretty(f->poly + i, vars, ctx);
		flint_printf(")^");
        fmpz_print(f->exp + i);
    }
}

void fmpz_mod_mpolyu_print_pretty(
    const fmpz_mod_mpolyu_t poly,
    const char ** x,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fmpz_mod_bpoly_print_pretty(
    const fmpz_mod_bpoly_t A,
    const char * xvar,
    const char * yvar,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_mod_poly_is_zero(A->coeffs + i, ctx))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, yvar, ctx);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_polyu3_print_pretty(
    const fmpz_mod_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        fmpz_print(A->coeffs + i);
        flint_printf("*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_mpolyv_print_pretty(
    const fmpz_mod_mpolyv_t poly,
    const char ** x,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        fmpz_mod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}

void fmpz_mod_mpolyn_print_pretty(
    const fmpz_mod_mpolyn_t A,
    const char ** x_in,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_poly_struct * coeff = A->coeffs;
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
        fmpz_mod_poly_print_pretty(coeff + i, "v", ctx->ffinfo);
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
