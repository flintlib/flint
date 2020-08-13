/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* evaluate B(xbar) at xbar = C */
int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A,
                     const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    slong i;
    fmpz_mat_t M;

    FLINT_ASSERT(A != B);

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctxAC);
        return 1;
    }

    fmpz_mat_init(M, ctxAC->minfo->nfields + 1, ctxB->minfo->nfields);
    fmpz_mat_zero(M);

    for (i = 0; i < ctxB->minfo->nvars; i++)
    {
        if (C[i]->length > 1)
            goto matrix_no_good;

        if (C[i]->length == 0)
        {
            mpoly_compose_mat_fill_column(M, NULL, 0, i,
                                                    ctxB->minfo, ctxAC->minfo);
        }
        else
        {
            if (!fmpz_is_one(C[i]->coeffs + 0))
                goto matrix_no_good;

            mpoly_compose_mat_fill_column(M, C[i]->exps, C[i]->bits, i,
                                                    ctxB->minfo, ctxAC->minfo);
        }
    }

    _fmpz_mpoly_compose_mat(A, B, M, ctxB, ctxAC);

    fmpz_mat_clear(M);

    return 1;

matrix_no_good:

    fmpz_mat_clear(M);

    for (i = 0; i < ctxB->minfo->nvars; i++)
    {
        if (C[i]->length > 1)
        {
            return fmpz_mpoly_compose_fmpz_mpoly_horner(A, B, C, ctxB, ctxAC);
        }
    }

    return fmpz_mpoly_compose_fmpz_mpoly_geobucket(A, B, C, ctxB, ctxAC);
}

