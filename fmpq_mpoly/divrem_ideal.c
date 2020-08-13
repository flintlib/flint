/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


/* Assumes divisor polys don't alias any output polys */
void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** Q, fmpq_mpoly_t R,
              const fmpq_mpoly_t A, fmpq_mpoly_struct * const * B, slong len,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t scale;
    fmpq_t t;
    fmpz_mpoly_struct ** Qarr, ** Barr;
    TMP_INIT;

    /* check none of the divisor polynomials is zero */
    for (i = 0; i < len; i++)
    {
        if (fmpq_mpoly_is_zero(B[i], ctx))
        {
            flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divrem_ideal");
        }
    }

    /* dividend is zero, write out quotients and remainder */
    if (fmpq_mpoly_is_zero(A, ctx))
    {
        for (i = 0; i < len; i++)
            fmpq_mpoly_zero(Q[i], ctx);
        fmpq_mpoly_zero(R, ctx);
        return;
    }

    TMP_START;
    Qarr = (fmpz_mpoly_struct **) TMP_ALLOC(len*sizeof(fmpz_mpoly_struct *));
    Barr = (fmpz_mpoly_struct **) TMP_ALLOC(len*sizeof(fmpz_mpoly_struct *));

    for (i = 0; i < len; i++)
    {
        Qarr[i] = Q[i]->zpoly;
        Barr[i] = B[i]->zpoly;
    }
    fmpz_init(scale);
    fmpz_mpoly_quasidivrem_ideal_heap(scale, Qarr, R->zpoly,
                                               A->zpoly, Barr, len, ctx->zctx);

    fmpq_init(t);
    fmpq_div_fmpz(t, A->content, scale);
    for (i = 0; i < len; i++)
        fmpq_div(Q[i]->content, t, B[i]->content);
    fmpq_swap(t, R->content);
    fmpq_clear(t);
    fmpz_clear(scale);

    for (i = 0; i < len; i++)
        fmpq_mpoly_reduce(Q[i], ctx);
    fmpq_mpoly_reduce(R, ctx);

    TMP_END;
}
