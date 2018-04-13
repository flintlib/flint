/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_pow_si(fmpq_mpoly_t qpoly1, const fmpq_mpoly_t qpoly2,
                                        slong pow, const fmpq_mpoly_ctx_t qctx)
{
    fmpz_mpoly_struct * poly1 = qpoly1->zpoly;
    const fmpz_mpoly_struct * poly2 = qpoly2->zpoly;
    const fmpz_mpoly_ctx_struct * ctx = qctx->zctx;

    if (pow < 0)
        flint_throw(FLINT_ERROR, "Negative power in fmpq_mpoly_pow_si");

    fmpq_pow_si(qpoly1->content, qpoly2->content, pow);
    fmpz_mpoly_pow_fps(poly1, poly2, pow, ctx);
    return;
}
