/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_div_fmpz(poly1->content, poly2->content, c);
    fmpz_mpoly_set(poly1->zpoly, poly2->zpoly, ctx->zctx);
}

void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, ulong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_ui(C, c);
    fmpq_mpoly_scalar_div_fmpz(poly1, poly2, C, ctx);
    fmpz_clear(C);
}

void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, slong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init(C);
    fmpz_set_si(C, c);
    fmpq_mpoly_scalar_div_fmpz(poly1, poly2, C, ctx);
    fmpz_clear(C);
}

