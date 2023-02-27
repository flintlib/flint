/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

/* evaluate B(x_1,...,x_n) at x_i = y_c[i], y_j are vars of ctxAC */
void fmpq_mpoly_compose_fmpq_mpoly_gen(fmpq_mpoly_t A,
                             const fmpq_mpoly_t B, const slong * c,
                     const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)
{
    fmpq_set(A->content, B->content);
    fmpz_mpoly_compose_fmpz_mpoly_gen(A->zpoly, B->zpoly, c,
                                                      ctxB->zctx, ctxAC->zctx);
    fmpq_mpoly_reduce(A, ctxAC);
    return;
}

