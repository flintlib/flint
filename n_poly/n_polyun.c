/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "mpn_extras.h"
#include "nmod_vec.h"


int n_polyun_mod_is_canonical(const n_polyun_t A, nmod_t mod)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_mod_is_canonical(A->terms[i].coeff, mod) ||
            n_poly_is_zero(A->terms[i].coeff))
        {
            return 0;
        }
        if (i > 0 && A->terms[i].exp >= A->terms[i - 1].exp)
            return 0;
    }
    return 1;
}

void n_polyun_clear(n_polyun_t A)
{
    slong i;

    for (i = 0; i < A->alloc; i++)
        n_poly_clear(A->terms[i].coeff);

    flint_free(A->terms);
}

void n_polyun_realloc(n_polyun_t A, slong len)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    A->terms = FLINT_ARRAY_REALLOC(A->terms, new_alloc, n_polyun_term_struct);

    for (i = old_alloc; i < new_alloc; i++)
        n_poly_init(A->terms[i].coeff);

    A->alloc = new_alloc;
}

void n_polyu1n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * varlast)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        n_poly_print_pretty(A->terms[i].coeff, varlast);
        flint_printf(")*%s^%wu",
            var0, A->terms[i].exp);
    }

    if (first)
        flint_printf("0");
}

void n_polyu2n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * varlast)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        n_poly_print_pretty(A->terms[i].coeff, varlast);
        flint_printf(")*%s^%wu*%s^%wu",
            var0, extract_exp(A->terms[i].exp, 1, 2),
            var1, extract_exp(A->terms[i].exp, 0, 2));
    }

    if (first)
        flint_printf("0");
}

void n_polyu3n_print_pretty(
    const n_polyun_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const char * varlast)
{
    slong i;
    int first = 1;

    for (i = 0; i < A->length; i++)
    {
        if (!first)
            printf(" + ");
        first = 0;
        flint_printf("(");
        n_poly_print_pretty(A->terms[i].coeff, varlast);
        flint_printf(")*%s^%wu*%s^%wu*%s^%wu",
            var0, extract_exp(A->terms[i].exp, 2, 3),
            var1, extract_exp(A->terms[i].exp, 1, 3),
            var2, extract_exp(A->terms[i].exp, 0, 3));
    }

    if (first)
        flint_printf("0");
}


void n_polyun_set(n_polyun_t A, const n_polyun_t B)
{
    slong i;
    n_polyun_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        A->terms[i].exp = B->terms[i].exp;
        n_poly_set(A->terms[i].coeff, B->terms[i].coeff);
    }
    A->length = B->length;
}

int n_polyun_equal(
    const n_polyun_t A,
    const n_polyun_t B)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (A->terms[i].exp != B->terms[i].exp)
            return 0;

        if (!n_poly_equal(A->terms[i].coeff, B->terms[i].coeff))
            return 0;
    }

    return 1;
}
