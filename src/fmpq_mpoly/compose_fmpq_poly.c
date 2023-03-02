/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_compose_fmpq_poly(fmpq_poly_t A,
                         const fmpq_mpoly_t B, fmpq_poly_struct * const * C,
                                                  const fmpq_mpoly_ctx_t ctxB)
{
    int success = 0;
    slong i;
    fmpq * scales;
    fmpz_poly_struct ** Czpoly;
    fmpz_mpoly_t newB;
    fmpq_t Acontent;
    fmpz_poly_t Azpoly;
    slong nvarsB = ctxB->zctx->minfo->nvars;
    TMP_INIT;

    if (B->zpoly->length == 0)
    {
        fmpq_poly_zero(A);
        return 1;
    }

    TMP_START;

    fmpq_init(Acontent);
    fmpz_poly_init(Azpoly);

    Czpoly = (fmpz_poly_struct **) TMP_ALLOC(nvarsB*sizeof(fmpz_poly_struct *));

    /*
        scale B by the contents of the polynomials in C
        We are only borrowing B to feed it to fmpz_mpoly_compose.
    */
    scales = (fmpq *) TMP_ALLOC(nvarsB*sizeof(fmpq));
    for (i = 0; i < nvarsB; i++)
    {
        /*
            since the fmpq_poly_t does not have a fmpz_poly_t in it,
            we have to manually copy the relevant struct members
        */
        Czpoly[i] = (fmpz_poly_struct *) flint_malloc(sizeof(fmpz_poly_struct));
        Czpoly[i]->coeffs = C[i]->coeffs;
        Czpoly[i]->alloc = C[i]->alloc;
        Czpoly[i]->length = C[i]->length;
        /* and manually set the scales */
        fmpz_init_set_ui(fmpq_numref(scales + i), UWORD(1));
        fmpz_init_set(fmpq_denref(scales + i), C[i]->den);
    }

    *newB = *B->zpoly;
    newB->coeffs = _fmpz_vec_init(B->zpoly->length);

    if (!_fmpq_mpoly_rescale(Acontent, newB->coeffs, B, scales, ctxB))
        goto cleanup;

    if (!fmpz_mpoly_compose_fmpz_poly(Azpoly, newB, Czpoly, ctxB->zctx))
        goto cleanup;

    fmpq_poly_set_fmpz_poly(A, Azpoly);
    fmpq_poly_scalar_mul_fmpq(A, A, Acontent);

    success = 1;

cleanup:

    _fmpz_vec_clear(newB->coeffs, B->zpoly->length);

    for (i = 0; i < nvarsB; i++)
    {
        flint_free(Czpoly[i]);
        fmpq_clear(scales + i);
    }

    fmpq_clear(Acontent);
    fmpz_poly_clear(Azpoly);

    TMP_END;

    return success;
}


