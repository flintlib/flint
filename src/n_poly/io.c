/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fq_nmod.h"
#include "n_poly.h"
#include "mpoly.h"

void n_polyu3_print_pretty(
    const n_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("%wu*%s^%wu*%s^%wu*%s^%wu",
            A->coeffs[i],
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void n_fq_bpoly_print_pretty(
    const n_fq_bpoly_t A,
    const char * xvar,
    const char * yvar,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (i + 1 != A->length && n_fq_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_fq_poly_print_pretty(A->coeffs + i, yvar, ctx);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void n_bpoly_print_pretty(
    const n_bpoly_t A,
    const char * xvar,
    const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (i < A->length - 1 && n_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

char * n_fq_get_str_pretty(
    const mp_limb_t * a,
    const fq_nmod_ctx_t ctx)
{
    char * s;
    fq_nmod_t A;
    fq_nmod_init(A, ctx);
    n_fq_get_fq_nmod(A, a, ctx);
    s = fq_nmod_get_str_pretty(A, ctx);
    fq_nmod_clear(A, ctx);
    return s;
}

int n_fq_fprint_pretty(
    FILE * file,
    const mp_limb_t * a,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    int first;

    first = 1;
    for (i = d - 1; i >= 0; i--)
    {
        if (a[i] == 0)
            continue;

        if (!first)
            flint_fprintf(file, "+");

        first = 0;
        flint_fprintf(file, "%wu", a[i]);

        if (i > 0)
        {
            flint_fprintf(file, "*%s", ctx->var);
            if (i > 1)
                flint_fprintf(file, "^%wd", i);
        }
    }

    if (first)
        flint_fprintf(file, "0");

    return 1;
}

void n_fq_print_pretty(const mp_limb_t * a, const fq_nmod_ctx_t ctx) { n_fq_fprint_pretty(stdout, a, ctx); }

void n_poly_print_pretty(const n_poly_t A, const char * x)
{
    slong i;
    int first = 1;

    for (i = A->length - 1; i >= 0; i--)
    {
        if (i < A->length - 1 && A->coeffs[i] == 0)
            continue;
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("%wu*%s^%wd", A->coeffs[i], x, i);
    }

    if (first)
        flint_printf("0");
}

void n_fq_poly_print_pretty(
    const n_fq_poly_t A,
    const char * x,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (i + 1 != A->length && _n_fq_is_zero(A->coeffs + d*i, d))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_fq_print_pretty(A->coeffs + d*i, ctx);
        flint_printf(")*%s^%wd", x, i);
    }

    if (first)
        flint_printf("0");
}

void n_polyu1n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * varlast)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        n_poly_print_pretty(A->coeffs + i, varlast);
        flint_printf(")*%s^%wu", var0, A->exps[i]);
    }

    if (first)
        flint_printf("0");
}

void n_polyu2n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        n_poly_print_pretty(A->coeffs + i, varlast);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 1, 2),
            var1, extract_exp(A->exps[i], 0, 2));
    }

    if (first)
        flint_printf("0");
}

void n_polyu3n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            flint_printf(" + ");
        first = 0;
        flint_printf("(");
        n_poly_print_pretty(A->coeffs + i, varlast);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}
