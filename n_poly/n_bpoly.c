/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_bpoly_clear(n_bpoly_t A)
{
    slong i;
    if (A->alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        for (i = 0; i < A->alloc; i++)
            n_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
    }
}


void n_bpoly_realloc(n_bpoly_t A, slong len)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);

    if (len <= A->alloc)
        return;

    if (A->alloc > 0)
    {
        A->coeffs = (n_poly_struct *) flint_realloc(A->coeffs,
                                            new_alloc * sizeof(n_poly_struct));
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        A->coeffs = (n_poly_struct *) flint_malloc(
                                            new_alloc * sizeof(n_poly_struct));
    }

    for (i = old_alloc; i < new_alloc; i++)
        n_poly_init(A->coeffs + i);

    A->alloc = len;
}


void n_bpoly_print_pretty(
    const n_bpoly_t A,
    const char * xvar,
    const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (i < A->length - 1 && n_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}


slong n_bpoly_degree1(const n_bpoly_t A)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;    
}

int n_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_equal(A->coeffs + i, B->coeffs + i))
            return 0;
    }

    return 1;
}


void _n_bpoly_set(n_bpoly_t A, const n_bpoly_t B)
{
    slong i;

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        n_poly_set(A->coeffs + i, B->coeffs + i);
}

void n_bpoly_set_coeff_nonzero(n_bpoly_t A, slong xi, slong yi, mp_limb_t c)
{
    slong i;

    FLINT_ASSERT(c != 0);

    if (xi >= A->length)
    {
        n_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            n_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    n_poly_set_coeff_nonzero(A->coeffs + xi, yi, c);
    FLINT_ASSERT(!n_poly_is_zero(A->coeffs + A->length - 1));
}

void n_bpoly_set_coeff(n_bpoly_t A, slong xi, slong yi, mp_limb_t c)
{
    slong i;

    if (xi >= A->length)
    {
        n_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            n_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    n_poly_set_coeff(A->coeffs + xi, yi, c);
    while (A->length > 0 && n_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

void n_bpoly_set_poly_gen1(n_bpoly_t A, const n_poly_t B)
{
    n_bpoly_fit_length(A, 1);
	n_poly_set(A->coeffs + 0, B);
	A->length = !n_poly_is_zero(A->coeffs + 0);
}


void n_bpoly_set_poly_gen0(n_bpoly_t A, const n_poly_t B)
{
    slong i;
    n_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        n_poly_set_ui(A->coeffs + i, B->coeffs[i]);
    A->length = B->length;
}


void n_bpoly_one(n_bpoly_t A)
{
    n_bpoly_fit_length(A, 1);
    A->length = 1;
    n_poly_one(A->coeffs + 0);
}

