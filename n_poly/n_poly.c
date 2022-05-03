/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "flint-impl.h"
#include "n_poly.h"

int n_poly_is_canonical(const n_poly_t A)
{
    if (A->length < 0)
        return 0;

    return A->length < 1 || A->coeffs[A->length - 1] != 0;
}

void n_poly_realloc(n_poly_t A, slong len)
{
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    if (old_alloc > 0)
    {
        FLINT_ASSERT(A->coeffs != NULL);
        A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs,
                                                  new_alloc*sizeof(mp_limb_t));
    }
    else
    {
        FLINT_ASSERT(A->coeffs == NULL);
        A->coeffs = (mp_limb_t *) flint_malloc(new_alloc*sizeof(mp_limb_t));
    }
    A->alloc = new_alloc;
}


void n_poly_print_pretty(const n_poly_t A, const char * x)
{
    slong i;
    int first = 1;

    for (i = A->length - 1; i >= 0; i--)
    {
        if (i < A->length - 1 && A->coeffs[i] == 0)
            continue;
        if (!first)
            printf(" + ");
        first = 0;
        printf(WORD_FMT "u*%s^" WORD_FMT "d", A->coeffs[i], x, i);
    }

    if (first)
        printf("0");
}

void n_poly_set_coeff(n_poly_t poly, slong j, ulong c)
{
    n_poly_fit_length(poly, j + 1);

    if (j + 1 < poly->length) /* interior */
    {
        poly->coeffs[j] = c;
    }
    else if (j + 1 == poly->length) /* leading coeff */
    {
        if (c != 0)
        {
            poly->coeffs[j] = c;
        }
        else
        {
            poly->length--;
            _n_poly_normalise(poly);
        }
    }
    else if (c != 0) /* extend polynomial */
    {
        FLINT_MPN_ZERO(poly->coeffs + poly->length, j - poly->length);
        poly->coeffs[j] = c;
        poly->length = j + 1;
    }
}

