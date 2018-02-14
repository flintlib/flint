/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_mul(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                          const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(poly2, ctx) || fmpq_mpoly_is_zero(poly3, ctx))
    {
        fmpq_mpoly_zero(poly1, ctx);
        return;
    }

    fmpq_mul(poly1->content, poly2->content, poly3->content);
    fmpz_mpoly_mul_johnson(poly1->zpoly, poly2->zpoly, poly3->zpoly, ctx->zctx);
}
