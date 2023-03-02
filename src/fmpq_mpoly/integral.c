/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_integral(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                         slong var, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t s;
    fmpz_init(s);
    fmpz_mpoly_integral(poly1->zpoly, s, poly2->zpoly, var, ctx->zctx);
    fmpq_div_fmpz(poly1->content, poly2->content, s);
    fmpq_mpoly_reduce(poly1, ctx);
    fmpz_clear(s);
}
