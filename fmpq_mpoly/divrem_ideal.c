/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


/* Assumes divisor polys don't alias any output polys */
void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** q, fmpq_mpoly_t r,
    const fmpq_mpoly_t a, fmpq_mpoly_struct * const * b, slong len,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t scale;
    fmpq_t t;
    fmpz_mpoly_struct ** qarr, ** barr;
    TMP_INIT;


    /* check none of the divisor polynomials is zero */
    for (i = 0; i < len; i++)
    {
        if (fmpq_mpoly_is_zero(b[i], ctx))
            flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divrem_ideal");
    }

    /* dividend is zero, write out quotients and remainder */
    if (fmpq_mpoly_is_zero(a, ctx))
    {
        for (i = 0; i < len; i++)
            fmpq_mpoly_zero(q[i], ctx);
        fmpq_mpoly_zero(r, ctx);
        return;
    }

    TMP_START;
    qarr = (fmpz_mpoly_struct **) TMP_ALLOC(len*sizeof(fmpz_mpoly_struct *));
    barr = (fmpz_mpoly_struct **) TMP_ALLOC(len*sizeof(fmpz_mpoly_struct *));

    for (i = 0; i < len; i++)
    {
        qarr[i] = q[i]->zpoly;
        barr[i] = b[i]->zpoly;
    }
    fmpz_init(scale);
    fmpz_mpoly_quasidivrem_ideal_heap(scale, qarr, r->zpoly,
                                               a->zpoly, barr, len, ctx->zctx);

    fmpq_init(t);
    fmpq_div_fmpz(t, a->content, scale);
    for (i = 0; i < len; i++)
        fmpq_div(q[i]->content, t, b[i]->content);
    fmpq_swap(t, r->content);
    fmpq_clear(t);
    fmpz_clear(scale);

    for (i = 0; i < len; i++)
        fmpq_mpoly_canonicalise(q[i], ctx);
    fmpq_mpoly_canonicalise(r, ctx);

    TMP_END;
}
