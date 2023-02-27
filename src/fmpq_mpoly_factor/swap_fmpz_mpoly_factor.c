/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_mpoly_factor.h"

/* f = g*c, clobber g */
void _fmpq_mpoly_factor_swap_fmpz_mpoly_factor(
    fmpq_mpoly_factor_t f,
    fmpz_mpoly_factor_t g,
    const fmpq_t c,
    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
	fmpq_mpoly_factor_fit_length(f, g->num, ctx);
	for (i = 0; i < g->num; i++)
	{
		fmpz_swap(f->exp + i, g->exp + i);
		fmpq_one(f->poly[i].content);
		fmpz_mpoly_swap(f->poly[i].zpoly, g->poly + i, ctx->zctx);
        fmpq_mpoly_reduce(f->poly + i, ctx); /* just in case */
	}
	f->num = g->num;
	fmpq_mul_fmpz(f->constant, c, g->constant);
}
