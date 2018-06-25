/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_set_fmpq(fmpq_mpoly_t poly,
                                    const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_is_zero(c))
    {
        fmpq_mpoly_zero(poly, ctx);
        return;
    }

    fmpq_set(poly->content, c);
    fmpz_mpoly_one(poly->zpoly, ctx->zctx);
}
