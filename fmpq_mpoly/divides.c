/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


/* return 1 if quotient is exact */
int fmpq_mpoly_divides(fmpq_mpoly_t poly1,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    int res;

    if (fmpq_mpoly_is_zero(poly3, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divides");
    }

    if (fmpq_mpoly_is_zero(poly2, ctx))
    {
        fmpq_mpoly_zero(poly1, ctx);
        return 1;
    }

    res = fmpz_mpoly_divides_monagan_pearce(poly1->zpoly, poly2->zpoly, poly3->zpoly, ctx->zctx);
    if (!res)
    {
        fmpq_mpoly_zero(poly1, ctx);
        return 0;
    }

    fmpq_div(poly1->content, poly2->content, poly3->content);
    return 1;
}
