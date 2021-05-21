/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_mpoly.h"

void nmod_poly_stack_init(nmod_poly_stack_t S, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    S->bits = bits;
    S->ctx = ctx;

    S->poly_alloc = 0;
    S->poly_array = NULL;
    S->poly_top = 0;

    S->mpolyun_alloc = 0;
    S->mpolyun_array = NULL;
    S->mpolyun_top = 0;

    S->mpolyn_alloc = 0;
    S->mpolyn_array = NULL;
    S->mpolyn_top = 0;
}

void nmod_poly_stack_clear(nmod_poly_stack_t S)
{
    slong i;

    FLINT_ASSERT(S->poly_top == 0);

    for (i = 0; i < S->poly_alloc; i++)
    {
        n_poly_clear(S->poly_array[i]);
        flint_free(S->poly_array[i]);
    }
    if (S->poly_array)
        flint_free(S->poly_array);

    for (i = 0; i < S->mpolyun_alloc; i++)
    {
        nmod_mpolyun_clear(S->mpolyun_array[i], S->ctx);
        flint_free(S->mpolyun_array[i]);
    }
    if (S->mpolyun_array)
        flint_free(S->mpolyun_array);

    for (i = 0; i < S->mpolyn_alloc; i++)
    {
        nmod_mpolyn_clear(S->mpolyn_array[i], S->ctx);
        flint_free(S->mpolyn_array[i]);
    }
    if (S->mpolyn_array)
        flint_free(S->mpolyn_array);

    S->ctx = NULL;
}

void nmod_poly_stack_set_ctx(nmod_poly_stack_t S, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(S->ctx->minfo->nvars == ctx->minfo->nvars);

    S->ctx = ctx;

    for (i = 0; i < S->mpolyun_alloc; i++)
    {
        nmod_mpolyun_set_mod(S->mpolyun_array[i], S->ctx->mod);
    }

    for (i = 0; i < S->mpolyn_alloc; i++)
    {
        nmod_mpolyn_set_mod(S->mpolyn_array[i], S->ctx->mod);
    }
}

/* insure that k slots are available after top and return pointer to top */
n_poly_struct ** nmod_poly_stack_fit_request_poly(nmod_poly_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->poly_alloc >= S->poly_top);

    if (S->poly_top + k > S->poly_alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->poly_top + k);

        if (S->poly_array)
        {
            S->poly_array = (n_poly_struct **) flint_realloc(S->poly_array,
                                           newalloc*sizeof(n_poly_struct*));
        }
        else
        {
            S->poly_array = (n_poly_struct **) flint_malloc(
                                           newalloc*sizeof(n_poly_struct*));
        }

        for (i = S->poly_alloc; i < newalloc; i++)
        {
            S->poly_array[i] = (n_poly_struct *) flint_malloc(
                                                     sizeof(n_poly_struct));
            n_poly_init(S->poly_array[i]);
        }
        S->poly_alloc = newalloc;
    }

    return S->poly_array + S->poly_top;
}



nmod_mpolyun_struct ** nmod_poly_stack_fit_request_mpolyun(nmod_poly_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->mpolyun_alloc >= S->mpolyun_top);

    if (S->mpolyun_top + k > S->mpolyun_alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->mpolyun_top + k);

        if (S->mpolyun_array)
        {
            S->mpolyun_array = (nmod_mpolyun_struct **) flint_realloc(S->mpolyun_array,
                                           newalloc*sizeof(nmod_mpolyun_struct*));
        }
        else
        {
            S->mpolyun_array = (nmod_mpolyun_struct **) flint_malloc(
                                           newalloc*sizeof(nmod_mpolyun_struct*));
        }

        for (i = S->mpolyun_alloc; i < newalloc; i++)
        {
            S->mpolyun_array[i] = (nmod_mpolyun_struct *) flint_malloc(
                                                     sizeof(nmod_mpolyun_struct));
            nmod_mpolyun_init(S->mpolyun_array[i], S->bits, S->ctx);
        }
        S->mpolyun_alloc = newalloc;
    }

    return S->mpolyun_array + S->mpolyun_top;
}

nmod_mpolyn_struct ** nmod_poly_stack_fit_request_mpolyn(nmod_poly_stack_t S, slong k)
{
    slong newalloc, i;

    FLINT_ASSERT(S->mpolyn_alloc >= S->mpolyn_top);

    if (S->mpolyn_top + k > S->mpolyn_alloc)
    {
        newalloc = FLINT_MAX(WORD(1), S->mpolyn_top + k);

        if (S->mpolyn_array)
        {
            S->mpolyn_array = (nmod_mpolyn_struct **) flint_realloc(S->mpolyn_array,
                                           newalloc*sizeof(nmod_mpolyn_struct*));
        }
        else
        {
            S->mpolyn_array = (nmod_mpolyn_struct **) flint_malloc(
                                           newalloc*sizeof(nmod_mpolyn_struct*));
        }

        for (i = S->mpolyn_alloc; i < newalloc; i++)
        {
            S->mpolyn_array[i] = (nmod_mpolyn_struct *) flint_malloc(
                                                     sizeof(nmod_mpolyn_struct));
            nmod_mpolyn_init(S->mpolyn_array[i], S->bits, S->ctx);
        }
        S->mpolyn_alloc = newalloc;
    }

    return S->mpolyn_array + S->mpolyn_top;
}

