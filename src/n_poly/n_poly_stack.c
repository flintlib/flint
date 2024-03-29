/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_poly_stack_init(n_poly_stack_t S)
{
    S->alloc = 0;
    S->array = NULL;
    S->top = 0;
}

void n_poly_stack_clear(n_poly_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        n_poly_clear(S->array[i]);
        flint_free(S->array[i]);
    }
    if (S->array)
        flint_free(S->array);
}

/* insure that k slots are available after top and return pointer to top */
n_poly_struct ** n_poly_stack_fit_request(n_poly_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->top + k);

        if (S->array)
        {
            S->array = (n_poly_struct **) flint_realloc(S->array,
                                              newalloc*sizeof(n_poly_struct*));
        }
        else
        {
            S->array = (n_poly_struct **) flint_malloc(
                                              newalloc*sizeof(n_poly_struct*));
        }

        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = (n_poly_struct *) flint_malloc(sizeof(n_poly_struct));
            n_poly_init(S->array[i]);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}
