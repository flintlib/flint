/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_content_vars(
    fmpq_mpoly_t g,
    const fmpq_mpoly_t A,
    slong * vars, slong num_vars,
    const fmpq_mpoly_ctx_t ctx)
{
    int success;

    success = fmpz_mpoly_content_vars(g->zpoly, A->zpoly, vars, num_vars, ctx->zctx);
    if (!success)
        return 0;

    if (g->zpoly->length > 0)
    {
        fmpz_one(fmpq_numref(g->content));
        fmpz_set(fmpq_denref(g->content), g->zpoly->coeffs + 0);
    }
    else
    {
        fmpq_zero(g->content);
    }

    return 1;
}

