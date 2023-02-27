/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* evaluate B(x_1,...,x_n) at x_i = y_c[i], y_j are vars of ctxAC */
void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const slong * c,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    fmpz_mat_t M;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctxAC);
        return;
    }

    fmpz_mat_init(M, ctxAC->minfo->nfields + 1, ctxB->minfo->nfields);
    mpoly_compose_mat_gen(M, c, ctxB->minfo, ctxAC->minfo);

    if (A == B)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init(T, ctxAC);
        _fmpz_mpoly_compose_mat(T, B, M, ctxB, ctxAC);
        fmpz_mpoly_swap(A, T, ctxAC);
        fmpz_mpoly_clear(T, ctxAC);
    }
    else
    {
        _fmpz_mpoly_compose_mat(A, B, M, ctxB, ctxAC);
    }

    fmpz_mat_clear(M);

    return;
}

