/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* evaluate B(x_1,...,x_n) at x_i = y_c[i], y_j are vars of ctxAC */
void fq_nmod_mpoly_compose_fq_nmod_mpoly_gen(fq_nmod_mpoly_t A,
                             const fq_nmod_mpoly_t B, const slong * c,
               const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)
{
    fmpz_mat_t M;

    if (B->length == 0)
    {
        fq_nmod_mpoly_zero(A, ctxAC);
        return;
    }

    fmpz_mat_init(M, ctxAC->minfo->nfields + 1, ctxB->minfo->nfields);
    mpoly_compose_mat_gen(M, c, ctxB->minfo, ctxAC->minfo);

    if (A == B)
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init(T, ctxAC);
        _fq_nmod_mpoly_compose_mat(T, B, M, ctxB, ctxAC);
        fq_nmod_mpoly_swap(A, T, ctxAC);
        fq_nmod_mpoly_clear(T, ctxAC);
    }
    else
    {
        _fq_nmod_mpoly_compose_mat(A, B, M, ctxB, ctxAC);
    }

    fmpz_mat_clear(M);

    return;
}

