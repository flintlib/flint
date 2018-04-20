/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_mul(poly1->content, poly2->content, c);
    if (fmpq_is_zero(poly1->content)) {
        fmpz_mpoly_zero(poly1->zpoly, ctx->zctx);
    } else {
        fmpz_mpoly_set(poly1->zpoly, poly2->zpoly, ctx->zctx);
    }
}
