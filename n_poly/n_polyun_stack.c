/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_polyun_stack_init(n_polyun_stack_t S)
{
    S->alloc = 0;
    S->array = NULL;
    S->top = 0;
}

void n_polyun_stack_clear(n_polyun_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        n_polyun_clear(S->array[i]);
        flint_free(S->array[i]);
    }

    flint_free(S->array);
}

n_polyun_struct ** n_polyun_stack_fit_request(n_polyun_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->top + k);
        S->array = FLINT_ARRAY_REALLOC(S->array, newalloc, n_polyun_struct *);
        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = FLINT_ARRAY_ALLOC(1, n_polyun_struct);
            n_polyun_init(S->array[i]);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}

