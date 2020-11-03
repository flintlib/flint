/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


int n_fq_bpoly_is_canonical(const n_fq_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A->length < 0)
        return 0;

    if (A->length > A->alloc)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_poly_is_canonical(A->coeffs + i, ctx))
            return 0;
        if (i + 1 == A->length && n_fq_poly_is_zero(A->coeffs + i))
            return 0;
    }

    return 1;
}


void n_fq_bpoly_print_pretty(
    const n_fq_bpoly_t A,
    const char * xvar,
    const char * yvar,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (i + 1 != A->length && n_fq_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_fq_poly_print_pretty(A->coeffs + i, yvar, ctx);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}


void n_fq_bpoly_one(n_fq_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    n_fq_bpoly_fit_length(A, 1);
    A->length = 1;
    n_fq_poly_one(A->coeffs + 0, ctx);
}

void n_fq_bpoly_set(
    n_fq_bpoly_t A,
    const n_fq_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A == B)
        return;

    n_fq_bpoly_fit_length(A, B->length);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
        n_fq_poly_set(A->coeffs + i, B->coeffs + i, ctx);
}

int n_fq_bpoly_equal(
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_poly_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
    }

    return 1;
}


void n_fq_bpoly_get_coeff_n_fq(
    mp_limb_t * c,
    const n_bpoly_t A,
    slong e0,
    slong e1,
    const fq_nmod_ctx_t ctx)
{
    if (e0 >= A->length)
        _n_fq_zero(c, fq_nmod_ctx_degree(ctx));
    else
        n_fq_poly_get_coeff_n_fq(c, A->coeffs + e0, e1, ctx);
}

void n_fq_bpoly_set_coeff_n_fq(
    n_bpoly_t A,
    slong e0,
    slong e1,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (e0 >= A->length)
    {
        n_fq_bpoly_fit_length(A, e0 + 1);
        for (i = A->length; i <= e0; i++)
            n_fq_poly_zero(A->coeffs + i);
        A->length = e0 + 1;
    }

    n_fq_poly_set_coeff_n_fq(A->coeffs + e0, e1, c, ctx);

    n_fq_bpoly_normalise(A);

    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx));
}


void n_fq_bpoly_get_coeff_fq_nmod(
    fq_nmod_t c,
    const n_bpoly_t A,
    slong e0,
    slong e1,
    const fq_nmod_ctx_t ctx)
{
    if (e0 >= A->length)
        fq_nmod_zero(c, ctx);
    else
        n_fq_poly_get_coeff_fq_nmod(c, A->coeffs + e0, e1, ctx);
}


void n_fq_bpoly_set_fq_nmod_poly_gen0(
    n_fq_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_bpoly_fit_length(A, B->length);
    A->length = 0;
    for (i = 0; i < B->length; i++)
        n_fq_poly_set_fq_nmod(A->coeffs + i, B->coeffs + i, ctx);
    A->length = B->length;
    n_bpoly_normalise(A);
}

void n_fq_bpoly_set_n_fq_poly_gen0(
    n_fq_bpoly_t A,
    const n_fq_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    n_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        n_fq_poly_set_n_fq(A->coeffs + i, B->coeffs + d*i, ctx);
    A->length = B->length;
    n_bpoly_normalise(A);
}

void n_fq_bpoly_set_n_fq_poly_gen1(
    n_fq_bpoly_t A,
    const n_fq_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    n_bpoly_fit_length(A, 1);
	n_fq_poly_set(A->coeffs + 0, B, ctx);
	A->length = !n_poly_is_zero(A->coeffs + 0);
}


void n_fq_bpoly_derivative_gen0(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;

    if (Blen < 2)
    {
        n_bpoly_zero(A);
        return;
    }

    n_bpoly_fit_length(A, Blen - 1);

    for (i = 1; i < Blen; i++)
        n_fq_poly_scalar_mul_ui(A->coeffs + i - 1, B->coeffs + i, i, ctx);

    A->length = Blen - 1;
    n_bpoly_normalise(A);
}


void n_fq_bpoly_scalar_mul_n_fq(
    n_fq_bpoly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    if (_n_fq_is_zero(c, d))
    {
        A->length = 0;
        return;
    }

    if (_n_fq_is_one(c, d))
    {
        return;
    }

    for (i = 0; i < A->length; i++)
        n_fq_poly_scalar_mul_n_fq(A->coeffs + i, A->coeffs + i, c, ctx);
}
