/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A,
                   const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C,
                     const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)
{
    int success = 0;
    slong i;
    fmpq * scales;
    fmpz_mpoly_struct ** Czpoly;
    fmpz_mpoly_t newB;
    slong nvarsB = ctxB->zctx->minfo->nvars;
    TMP_INIT;

    if (fmpq_mpoly_is_zero(B, ctxB))
    {
        fmpq_mpoly_zero(A, ctxAC);
        return 1;
    }

    TMP_START;

    Czpoly = (fmpz_mpoly_struct **) TMP_ALLOC(nvarsB*sizeof(fmpz_mpoly_struct *));

    /*
        scale B by the contents of the polynomials in C
        We are only borrowing B to feed it to fmpz_mpoly_compose.
        There might be zero coeffs in newBcoeffs, but fmpz_mpoly_compose
            should have no problem with zero coeffs.
    */
    scales = (fmpq *) TMP_ALLOC(nvarsB*sizeof(fmpq));
    for (i = 0; i < nvarsB; i++)
    {
        Czpoly[i] = C[i]->zpoly;
        /* we are only borrowing the content of each of the C[i] */
        *(scales + i) = *C[i]->content;
    }
    *newB = *B->zpoly;
    newB->coeffs = _fmpz_vec_init(B->zpoly->length);

    if (!_fmpq_mpoly_rescale(A->content, newB->coeffs, B, scales, ctxB))
        goto cleanup;

    if (!fmpz_mpoly_compose_fmpz_mpoly(A->zpoly, newB, Czpoly,
                                                      ctxB->zctx, ctxAC->zctx))
        goto cleanup;

    fmpq_mpoly_reduce(A, ctxAC);

    success = 1;

cleanup:

    TMP_END;

    _fmpz_vec_clear(newB->coeffs, B->zpoly->length);

    return success;
}

