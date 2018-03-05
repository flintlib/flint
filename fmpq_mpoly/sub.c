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
void fmpq_mpoly_sub(fmpq_mpoly_t poly1,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t temp2, temp3;
    fmpz_t n2d3, d2n3, d2d3, one;

    fmpz_init(n2d3);
    fmpz_init(d2n3);
    fmpz_init(d2d3);
    fmpz_init_set_ui(one, 1);
    fmpz_mpoly_init(temp2, ctx->zctx);
    fmpz_mpoly_init(temp3, ctx->zctx);

    fmpz_mul(n2d3, fmpq_numref(poly2->content), fmpq_denref(poly3->content));
    fmpz_mul(d2n3, fmpq_denref(poly2->content), fmpq_numref(poly3->content));
    fmpz_mul(d2d3, fmpq_denref(poly2->content), fmpq_denref(poly3->content));

    fmpz_mpoly_scalar_mul_fmpz(temp2, poly2->zpoly, n2d3, ctx->zctx);
    fmpz_mpoly_scalar_mul_fmpz(temp3, poly3->zpoly, d2n3, ctx->zctx);

    fmpz_mpoly_sub(poly1->zpoly, temp2, temp3, ctx->zctx);
    fmpq_set_fmpz_frac(poly1->content, one, d2d3);

    fmpz_mpoly_clear(temp3, ctx->zctx);
    fmpz_mpoly_clear(temp2, ctx->zctx);
    fmpz_clear(one);
    fmpz_clear(d2d3);
    fmpz_clear(d2n3);
    fmpz_clear(n2d3);

    fmpq_mpoly_canonicalise(poly1, ctx);
}
