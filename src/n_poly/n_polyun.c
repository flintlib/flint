/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"

int n_polyun_mod_is_canonical(const n_polyun_t A, nmod_t mod)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_mod_is_canonical(A->coeffs + i, mod) ||
            n_poly_is_zero(A->coeffs + i))
        {
            return 0;
        }
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}

void n_polyun_clear(n_polyun_t A)
{
    slong i;

    for (i = 0; i < A->alloc; i++)
        n_poly_clear(A->coeffs + i);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void n_polyun_realloc(n_polyun_t A, slong len)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);
    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, n_poly_struct);

    for (i = old_alloc; i < new_alloc; i++)
        n_poly_init(A->coeffs + i);

    A->alloc = new_alloc;
}

void n_polyun_set(n_polyun_t A, const n_polyun_t B)
{
    slong i;
    n_polyun_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        n_poly_set(A->coeffs + i, B->coeffs + i);
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
        if (A->exps[i] != B->exps[i])
            return 0;

        if (!n_poly_equal(A->coeffs + i, B->coeffs + i))
            return 0;
    }

    return 1;
}
