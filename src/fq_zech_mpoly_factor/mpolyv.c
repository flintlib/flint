/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


void fq_zech_mpolyv_clear(fq_zech_mpolyv_t A, const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    if (A->alloc > 0)
    {
        for (i = 0; i < A->alloc; i++)
            fq_zech_mpoly_clear(A->coeffs + i, ctx);
        flint_free(A->coeffs);
    }
}

void fq_zech_mpolyv_fit_length(
    fq_zech_mpolyv_t A,
    slong length,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    if (old_alloc > 0)
    {
        A->coeffs = (fq_zech_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_zech_mpoly_struct));
    }
    else
    {
        A->coeffs = (fq_zech_mpoly_struct *) flint_malloc(
                                       new_alloc*sizeof(fq_zech_mpoly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        fq_zech_mpoly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}

void fq_zech_mpolyv_set_coeff(
    fq_zech_mpolyv_t A,
    slong i,
    fq_zech_mpoly_t c, /* clobbered */
    const fq_zech_mpoly_ctx_t ctx)
{
    slong j;
    FLINT_ASSERT(!fq_zech_mpoly_is_zero(c, ctx));
    fq_zech_mpolyv_fit_length(A, i + 1, ctx);
    for (j = A->length; j < i; j++)
        fq_zech_mpoly_zero(A->coeffs + j, ctx);
    fq_zech_mpoly_swap(A->coeffs + i, c, ctx);
    A->length = FLINT_MAX(A->length, i + 1);
}


void fq_zech_mpoly_to_mpolyv(
    fq_zech_mpolyv_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_t xalpha,
    const fq_zech_mpoly_ctx_t ctx)
{
    fq_zech_mpoly_t Q, T;

    fq_zech_mpoly_init(Q, ctx);
    fq_zech_mpoly_init(T, ctx);

    fq_zech_mpolyv_fit_length(A, 8, ctx);
    fq_zech_mpoly_divrem(Q, A->coeffs + 0, B, xalpha, ctx);
    A->length = 1;

    while (!fq_zech_mpoly_is_zero(Q, ctx))
    {
        fq_zech_mpolyv_fit_length(A, A->length + 1, ctx);
        fq_zech_mpoly_divrem(T, A->coeffs + A->length, Q, xalpha, ctx);
        fq_zech_mpoly_swap(Q, T, ctx);
        A->length++;
    }

    while (A->length > 0 && fq_zech_mpoly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;

    fq_zech_mpoly_clear(Q, ctx);
    fq_zech_mpoly_clear(T, ctx);
}


void fq_zech_mpoly_from_mpolyv(
    fq_zech_mpoly_t A,
    const fq_zech_mpolyv_t B,
    const fq_zech_mpoly_t xalpha,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    fq_zech_mpoly_t T;

    fq_zech_mpoly_init(T, ctx);

    fq_zech_mpoly_zero(A, ctx);
    for (i = B->length - 1; i >= 0; i--)
    {
        fq_zech_mpoly_mul(T, A, xalpha, ctx);
        fq_zech_mpoly_add(A, T, B->coeffs + i, ctx);
    }

    fq_zech_mpoly_clear(T, ctx);
}
