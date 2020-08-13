/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c, const fmpq_mpoly_t poly1,
                         const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_one(fmpq_denref(c));
    fmpz_mpoly_get_coeff_fmpz_monomial(fmpq_numref(c), poly1->zpoly,
                                                      poly2->zpoly, ctx->zctx);
    fmpq_mul(c, c, poly1->content);
}
