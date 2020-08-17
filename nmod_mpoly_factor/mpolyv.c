/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


void nmod_mpolyv_clear(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
}


void nmod_mpolyv_print_pretty(
    const nmod_mpolyv_t poly,
    const char ** x,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < poly->length; i++)
    {
        flint_printf("coeff[%wd]: ", i);
        nmod_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf("\n");
    }
}


void nmod_mpolyv_fit_length(
    nmod_mpolyv_t A,
    slong length,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    if (old_alloc > 0)
    {
        A->coeffs = (nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpoly_struct));
    }
    else
    {
        A->coeffs = (nmod_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(nmod_mpoly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        nmod_mpoly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}

void nmod_mpolyv_set_coeff(
    nmod_mpolyv_t A,
    slong i,
    nmod_mpoly_t c, /* clobbered */
    const nmod_mpoly_ctx_t ctx)
{
    slong j;
    FLINT_ASSERT(!nmod_mpoly_is_zero(c, ctx));
    nmod_mpolyv_fit_length(A, i + 1, ctx);
    for (j = A->length; j < i; j++)
        nmod_mpoly_zero(A->coeffs + j, ctx);
    nmod_mpoly_swap(A->coeffs + i, c, ctx);
    A->length = FLINT_MAX(A->length, i + 1);
}


void nmod_mpoly_to_mpolyv(
    nmod_mpolyv_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t xalpha,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t Q, T;

    nmod_mpoly_init(Q, ctx);
    nmod_mpoly_init(T, ctx);

    nmod_mpolyv_fit_length(A, 8, ctx);
    nmod_mpoly_divrem(Q, A->coeffs + 0, B, xalpha, ctx);
    A->length = 1;

    while (!nmod_mpoly_is_zero(Q, ctx))
    {
        nmod_mpolyv_fit_length(A, A->length + 1, ctx);
        nmod_mpoly_divrem(T, A->coeffs + A->length, Q, xalpha, ctx);
        nmod_mpoly_swap(Q, T, ctx);
        A->length++;
    }

    while (A->length > 0 && nmod_mpoly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;

    nmod_mpoly_clear(Q, ctx);
    nmod_mpoly_clear(T, ctx);
}


void nmod_mpoly_from_mpolyv(
    nmod_mpoly_t A,
    const nmod_mpolyv_t B,
    const nmod_mpoly_t xalpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_mpoly_t T;

    nmod_mpoly_init(T, ctx);

    nmod_mpoly_zero(A, ctx);
    for (i = B->length - 1; i >= 0; i--)
    {
        nmod_mpoly_mul(T, A, xalpha, ctx);
        nmod_mpoly_add(A, T, B->coeffs + i, ctx);
    }

    nmod_mpoly_clear(T, ctx);
}


int _nmod_mpoly_vec_content_mpoly(
    nmod_mpoly_t g,
    const nmod_mpoly_struct * A,
    slong Alen,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpoly_zero(g, ctx);

    for (i = 0; i < Alen; i++)
    {
		if (!nmod_mpoly_gcd(g, g, A + i, ctx))
			return 0;
    }

    return 1;
}

