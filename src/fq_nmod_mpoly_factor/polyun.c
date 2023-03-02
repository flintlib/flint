/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


void n_polyu2n_fq_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        n_fq_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 1, 2),
            var1, extract_exp(A->exps[i], 0, 2));
    }

    if (first)
        flint_printf("0");
}

void n_polyu3n_fq_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        n_fq_poly_print_pretty(A->coeffs + i, varlast, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

int n_polyun_fq_is_canonical(
    const n_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_poly_is_canonical(A->coeffs + i, ctx))
            return 0;
        if (n_poly_is_zero(A->coeffs + i))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}
