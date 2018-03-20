/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void

fmpq_mpoly_sub_fmpq(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                    const fmpq_t x, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, fmpq_numref(poly2->content), fmpq_denref(x));
    fmpz_mpoly_scalar_mul_fmpz(poly1->zpoly, poly2->zpoly, t, ctx->zctx);
    fmpz_mul(t, fmpq_denref(poly2->content), fmpq_numref(x));
    fmpz_set_ui(fmpq_numref(poly1->content), 1);
    fmpz_mul(fmpq_denref(poly1->content), fmpq_denref(poly2->content), fmpq_denref(x));
    fmpz_mpoly_sub_fmpz(poly1->zpoly, poly1->zpoly, t, ctx->zctx);
    fmpz_clear(t);

    fmpq_mpoly_canonicalise(poly1, ctx);
    return;
}
