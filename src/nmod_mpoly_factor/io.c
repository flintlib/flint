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
#include "mpoly.h"
#include "nmod_mpoly_factor.h"

void nmod_mpoly_factor_print_pretty(const nmod_mpoly_factor_t f,
                                const char ** vars, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    flint_printf("%wu", f->constant);
    for (i = 0; i < f->num; i++)
    {
        flint_printf("\n*(", i);
        nmod_mpoly_print_pretty(f->poly + i, vars, ctx);
		flint_printf(")^");
        fmpz_print(f->exp + i);
    }
}

void nmod_mpolyu3_print_pretty(
    const nmod_mpolyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char ** vars,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        nmod_mpoly_print_pretty(A->coeffs + i, vars, ctx);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->exps[i], 2, 3),
            var1, extract_exp(A->exps[i], 1, 3),
            var2, extract_exp(A->exps[i], 0, 3));
    }

    if (first)
        flint_printf("0");
}

void nmod_mpolyv_print_pretty(
    const nmod_mpolyv_t poly,
    const char ** x,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}
