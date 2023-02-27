/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"


int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_factor_t zf;

	fmpz_mpoly_factor_init(zf, ctx->zctx);
	success = fmpz_mpoly_factor_squarefree(zf, A->zpoly, ctx->zctx);
    _fmpq_mpoly_factor_swap_fmpz_mpoly_factor(f, zf, A->content, ctx);
	fmpz_mpoly_factor_clear(zf, ctx->zctx);

	return success;
}

