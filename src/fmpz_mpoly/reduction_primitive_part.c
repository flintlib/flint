/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void
fmpz_mpoly_reduction_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_vec_t I, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t scale;
    fmpz_mpoly_struct ** Q, ** B;
    slong i, len;

    len = I->length;

    fmpz_init(scale);
    Q = flint_malloc(sizeof(fmpz_mpoly_struct *) * len);
    B = flint_malloc(sizeof(fmpz_mpoly_struct *) * len);

    for (i = 0; i < len; i++)
    {
        Q[i] = flint_malloc(sizeof(fmpz_mpoly_struct));
        fmpz_mpoly_init(Q[i], ctx);
        B[i] = fmpz_mpoly_vec_entry(I, i);
    }

    fmpz_mpoly_quasidivrem_ideal(scale, Q, res, f, B, len, ctx);
    fmpz_mpoly_primitive_part(res, res, ctx);

    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_clear(Q[i], ctx);
        flint_free(Q[i]);
    }

    flint_free(Q);
    flint_free(B);
    fmpz_clear(scale);
}
