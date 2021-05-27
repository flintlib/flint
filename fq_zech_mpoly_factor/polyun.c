/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


void fq_zech_polyun_clear(fq_zech_polyun_t A, const fq_zech_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_zech_poly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_zech_polyun_realloc(fq_zech_polyun_t A, slong len, const fq_zech_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, fq_zech_poly_struct);
    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);

    for (i = old_alloc; i < new_alloc; i++)
        fq_zech_poly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
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

int fq_zech_polyun_is_canonical(
    const fq_zech_polyun_t A,
    const fq_zech_ctx_t ctx)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (fq_zech_poly_is_zero(A->coeffs + i, ctx))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}
