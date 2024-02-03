/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fq_zech.h"
#include "fq_zech_poly.h"
#include "mpoly.h"
#include "fq_zech_mpoly_factor.h"

void fq_zech_mpoly_factor_print_pretty(
    const fq_zech_mpoly_factor_t f,
    const char ** vars,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;

    flint_printf("(");
    fq_zech_print_pretty(f->constant, ctx->fqctx);
    flint_printf(")");
    for (i = 0; i < f->num; i++)
    {
        flint_printf("\n*(", i);
        fq_zech_mpoly_print_pretty(f->poly + i, vars, ctx);
		flint_printf(")^");
        fmpz_print(f->exp + i);
    }
}

void fq_zech_polyu3_print_pretty(
    const fq_zech_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fq_zech_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        fq_zech_print_pretty(A->coeffs + i, ctx);
        flint_printf("%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void fq_zech_bpoly_print_pretty(
    const fq_zech_bpoly_t A,
    const char * var0,
    const char * var1,
    const fq_zech_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fq_zech_poly_is_zero(A->coeffs + i, ctx))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fq_zech_poly_print_pretty(A->coeffs + i, var1, ctx);
        flint_printf(")*%s^%wd", var0, i);
    }

    if (first)
        flint_printf("0");
}

void fq_zech_mpolyv_print_pretty(
    const fq_zech_mpolyv_t poly,
    const char ** x,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        fq_zech_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}

void fq_zech_polyu2n_print_pretty(
    const fq_zech_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fq_zech_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fq_zech_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 1, 2),
            var1, extract_exp(A->exps[i], 0, 2));
    }

    if (first)
        flint_printf("0");
}

void fq_zech_polyu3n_print_pretty(
    const fq_zech_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fq_zech_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        fq_zech_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}
