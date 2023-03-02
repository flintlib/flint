/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

void fq_nmod_mpoly_factor_print_pretty(
    const fq_nmod_mpoly_factor_t f,
    const char ** vars,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    flint_printf("(");
    fq_nmod_print_pretty(f->constant, ctx->fqctx);
    flint_printf(")");
    for (i = 0; i < f->num; i++)
    {
        flint_printf("\n*(", i);
        fq_nmod_mpoly_print_pretty(f->poly + i, vars, ctx);
		flint_printf(")^");
        fmpz_print(f->exp + i);
    }
}
