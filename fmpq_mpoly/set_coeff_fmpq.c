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


fmpq_mpoly_set_coeff_fmpq(fmpq_mpoly_t poly,
                           slong n, const fmpq_t x, const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_is_zero(x))
    {
        fmpz_mpoly_set_coeff_fmpz(poly->zpoly, n, fmpq_numref(x), ctx->zctx);

    } else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mul(t, fmpq_numref(poly->content), fmpq_denref(x));
        fmpz_mpoly_scalar_mul_fmpz(poly->zpoly, poly->zpoly, t, ctx->zctx);
        fmpz_mul(t, fmpq_denref(poly->content), fmpq_numref(x));
        fmpz_set_ui(fmpq_numref(poly->content), 1);
        fmpz_mul(fmpq_denref(poly->content), fmpq_denref(poly->content), fmpq_denref(x));
        fmpz_mpoly_set_coeff_fmpz(poly->zpoly, n, t, ctx->zctx);
        fmpz_clear(t);
    }

    fmpq_mpoly_canonicalise(poly, ctx);
    return;
}
