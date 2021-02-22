/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


void fmpz_mod_poly_stack_init(fmpz_mod_poly_stack_t S)
{
    S->alloc = 0;
    S->array = NULL;
    S->top = 0;
}

void fmpz_mod_poly_stack_clear(fmpz_mod_poly_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        fmpz_mod_poly_clear(S->array[i], NULL);
        flint_free(S->array[i]);
    }

    flint_free(S->array);
}

/* insure that k slots are available after top and return pointer to top */
fmpz_mod_poly_struct ** fmpz_mod_poly_stack_fit_request(
    fmpz_mod_poly_stack_t S,
    slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = FLINT_MAX(1, S->top + k);

        S->array = FLINT_ARRAY_REALLOC(S->array, newalloc, fmpz_mod_poly_struct *);

        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = FLINT_ARRAY_ALLOC(1, fmpz_mod_poly_struct);
            fmpz_mod_poly_init(S->array[i], NULL);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}



void fmpz_mod_bpoly_stack_init(fmpz_mod_bpoly_stack_t S)
{
    S->alloc = 0;
    S->array = NULL;
    S->top = 0;
}

void fmpz_mod_bpoly_stack_clear(fmpz_mod_bpoly_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        fmpz_mod_bpoly_clear(S->array[i], NULL);
        flint_free(S->array[i]);
    }

    flint_free(S->array);
}

/* insure that k slots are available after top and return pointer to top */
fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_fit_request(
    fmpz_mod_bpoly_stack_t S,
    slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = FLINT_MAX(1, S->top + k);

        S->array = FLINT_ARRAY_REALLOC(S->array, newalloc, fmpz_mod_bpoly_struct *);

        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = FLINT_ARRAY_ALLOC(1, fmpz_mod_bpoly_struct);
            fmpz_mod_bpoly_init(S->array[i], NULL);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}


void fmpz_mod_polyun_stack_init(fmpz_mod_polyun_stack_t S)
{
    S->alloc = 0;
    S->array = NULL;
    S->top = 0;
}

void fmpz_mod_polyun_stack_clear(fmpz_mod_polyun_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        fmpz_mod_polyun_clear(S->array[i], NULL);
        flint_free(S->array[i]);
    }

    flint_free(S->array);
}

/* insure that k slots are available after top and return pointer to top */
fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_fit_request(
    fmpz_mod_polyun_stack_t S,
    slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = FLINT_MAX(1, S->top + k);

        S->array = FLINT_ARRAY_REALLOC(S->array, newalloc, fmpz_mod_polyun_struct *);

        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = FLINT_ARRAY_ALLOC(1, fmpz_mod_polyun_struct);
            fmpz_mod_polyun_init(S->array[i], NULL);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}





void fmpz_mod_mpolyn_stack_init(fmpz_mod_mpolyn_stack_t S,
                           flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)
{
    S->alloc = 0;
    S->array = NULL;
    S->top = 0;
    S->bits = bits;
}

void fmpz_mod_mpolyn_stack_clear(
    fmpz_mod_mpolyn_stack_t S,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(S->top == 0);

    for (i = 0; i < S->alloc; i++)
    {
        fmpz_mod_mpolyn_clear(S->array[i], ctx);
        flint_free(S->array[i]);
    }

    flint_free(S->array);
}

/* insure that k slots are available after top and return pointer to top */
fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_fit_request(
            fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong newalloc, i;

    FLINT_ASSERT(S->alloc >= S->top);

    if (S->top + k > S->alloc)
    {
        newalloc = FLINT_MAX(1, S->top + k);

        S->array = FLINT_ARRAY_REALLOC(S->array, newalloc, fmpz_mod_mpolyn_struct *);

        for (i = S->alloc; i < newalloc; i++)
        {
            S->array[i] = FLINT_ARRAY_ALLOC(1, fmpz_mod_mpolyn_struct);
            fmpz_mod_mpolyn_init(S->array[i], S->bits, ctx);
        }
        S->alloc = newalloc;
    }

    return S->array + S->top;
}
