/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_mpolyv_clear(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mpoly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
}


void fmpz_mpolyv_print_pretty(
    const fmpz_mpolyv_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}


void fmpz_mpolyv_fit_length(
    fmpz_mpolyv_t A,
    slong length,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    if (old_alloc > 0)
    {
        A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
    }
    else
    {
        A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        fmpz_mpoly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}


void fmpz_mpolyv_set_coeff(
    fmpz_mpolyv_t A,
    slong i,
    fmpz_mpoly_t c, /* clobbered */
    const fmpz_mpoly_ctx_t ctx)
{
    slong j;
    FLINT_ASSERT(!fmpz_mpoly_is_zero(c, ctx));
    fmpz_mpolyv_fit_length(A, i + 1, ctx);
    for (j = A->length; j < i; j++)
        fmpz_mpoly_zero(A->coeffs + j, ctx);
    fmpz_mpoly_swap(A->coeffs + i, c, ctx);
    A->length = FLINT_MAX(A->length, i + 1);
}


void fmpz_mpoly_to_mpolyv(
    fmpz_mpolyv_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t xalpha,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t Q, T;

    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(T, ctx);

    fmpz_mpolyv_fit_length(A, 8, ctx);
    fmpz_mpoly_divrem(Q, A->coeffs + 0, B, xalpha, ctx);
    A->length = 1;

    while (!fmpz_mpoly_is_zero(Q, ctx))
    {
        fmpz_mpolyv_fit_length(A, A->length + 1, ctx);
        fmpz_mpoly_divrem(T, A->coeffs + A->length, Q, xalpha, ctx);
        fmpz_mpoly_swap(Q, T, ctx);
        A->length++;
    }

    while (A->length > 0 && fmpz_mpoly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;

    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(T, ctx);
}


void fmpz_mpoly_from_mpolyv(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpolyv_t B,
    const fmpz_mpoly_t xalpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t T;

    fmpz_mpoly_init(T, ctx);

    fmpz_mpoly_zero(A, ctx);
    for (i = B->length - 1; i >= 0; i--)
    {
        fmpz_mpoly_mul(T, A, xalpha, ctx);
        fmpz_mpoly_add(A, T, B->coeffs + i, ctx);
    }

    fmpz_mpoly_clear(T, ctx);

    fmpz_mpoly_repack_bits_inplace(A, Abits, ctx);
}


int _fmpz_mpoly_vec_content_mpoly(
    fmpz_mpoly_t g,
    const fmpz_mpoly_struct * A,
    slong Alen,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j1, j2;

    if (Alen <= 1)
    {
        if (Alen == 1 && A[0].length > 0)
        {
            if (fmpz_sgn(A[0].coeffs + 0) < 0)
                fmpz_mpoly_neg(g, A + 0, ctx);
            else
                fmpz_mpoly_set(g, A + 0, ctx);
        }
        else
        {
            fmpz_mpoly_zero(g, ctx);
        }
        return 1;
    }

    j1 = 0;
    j2 = 1;
    for (i = 2; i < Alen; i++)
    {
        if (A[i].length < A[j1].length)
            j1 = i;
        else if (A[i].length < A[j2].length)
            j2 = i;
    }

    FLINT_ASSERT(j1 != j2);

    if (!fmpz_mpoly_gcd(g, A + j1, A + j2, ctx))
        return 0;

    for (i = 0; i < Alen; i++)
    {
        if (i == j1 || i == j2)
            continue;

        if (!fmpz_mpoly_gcd(g, g, A + i, ctx))
            return 0;
    }

    return 1;
}


void _fmpz_mpoly_vec_divexact_mpoly(
    fmpz_mpoly_struct * A, slong Alen,
    const fmpz_mpoly_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < Alen; i++)
        fmpz_mpoly_divexact(A + i, A + i, c, ctx);
}

void _fmpz_mpoly_vec_mul_mpoly(
    fmpz_mpoly_struct * A, slong Alen,
    const fmpz_mpoly_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < Alen; i++)
        fmpz_mpoly_mul(A + i, A + i, c, ctx);
}

