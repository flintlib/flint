/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_add_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                    const fmpz_t x, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mpoly_scalar_mul_fmpz(poly1->zpoly, poly2->zpoly,
                                       fmpq_numref(poly2->content), ctx->zctx);
    fmpz_one(fmpq_numref(poly1->content));
    fmpz_mul(t, fmpq_denref(poly2->content), x);
    fmpz_set(fmpq_denref(poly1->content), fmpq_denref(poly2->content));
    fmpz_mpoly_add_fmpz(poly1->zpoly, poly1->zpoly, t, ctx->zctx);

    fmpz_clear(t);
    fmpq_mpoly_canonicalise(poly1, ctx);
    return;
}

void fmpq_mpoly_add_ui(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                           ulong x, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t X;
    fmpz_init_set_ui(X, x);
    fmpq_mpoly_add_fmpz(poly1, poly2, X, ctx);
    fmpz_clear(X);
}

void fmpq_mpoly_add_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                           slong x, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t X;
    fmpz_init(X);
    fmpz_set_si(X, x);
    fmpq_mpoly_add_fmpz(poly1, poly2, X, ctx);
    fmpz_clear(X);
}
