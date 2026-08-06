/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void
fmpz_mod_mpoly_vec_init(fmpz_mod_mpoly_vec_t vec, slong len, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (len == 0)
    {
        vec->p = NULL;
        vec->length = 0;
        vec->alloc = 0;
    }
    else
    {
        slong i;
        vec->p = flint_malloc(sizeof(fmpz_mod_mpoly_struct) * len);
        for (i = 0; i < len; i++)
            fmpz_mod_mpoly_init(vec->p + i, ctx);
        vec->length = vec->alloc = len;
    }
}

void
fmpz_mod_mpoly_vec_print(const fmpz_mod_mpoly_vec_t F, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    flint_printf("[");
    for (i = 0; i < F->length; i++)
    {
        fmpz_mod_mpoly_print_pretty(F->p + i, NULL, ctx);
        if (i < F->length - 1)
            flint_printf(", ");
    }
    flint_printf("]");
}

void
fmpz_mod_mpoly_vec_swap(fmpz_mod_mpoly_vec_t x, fmpz_mod_mpoly_vec_t y, const fmpz_mod_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    fmpz_mod_mpoly_vec_t tmp;
    *tmp = *x;
    *x = *y;
    *y = *tmp;
}

void
fmpz_mod_mpoly_vec_fit_length(fmpz_mod_mpoly_vec_t vec, slong len, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (len > vec->alloc)
    {
        slong i;

        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->p = flint_realloc(vec->p, len * sizeof(fmpz_mod_mpoly_struct));

        for (i = vec->alloc; i < len; i++)
            fmpz_mod_mpoly_init(vec->p + i, ctx);

        vec->alloc = len;
    }
}

void
fmpz_mod_mpoly_vec_clear(fmpz_mod_mpoly_vec_t vec, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < vec->alloc; i++)
        fmpz_mod_mpoly_clear(vec->p + i, ctx);

    flint_free(vec->p);
}

void
fmpz_mod_mpoly_vec_set(fmpz_mod_mpoly_vec_t dest, const fmpz_mod_mpoly_vec_t src, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (dest != src)
    {
        slong i;

        fmpz_mod_mpoly_vec_fit_length(dest, src->length, ctx);

        for (i = 0; i < src->length; i++)
            fmpz_mod_mpoly_set(dest->p + i, src->p + i, ctx);

        dest->length = src->length;
    }
}

void
fmpz_mod_mpoly_vec_append(fmpz_mod_mpoly_vec_t vec, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_vec_fit_length(vec, vec->length + 1, ctx);
    fmpz_mod_mpoly_set(vec->p + vec->length, f, ctx);
    vec->length++;
}

slong
fmpz_mod_mpoly_vec_insert_unique(fmpz_mod_mpoly_vec_t vec, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < vec->length; i++)
    {
        if (fmpz_mod_mpoly_equal(vec->p + i, f, ctx))
            return i;
    }

    fmpz_mod_mpoly_vec_append(vec, f, ctx);
    return vec->length - 1;
}

void
fmpz_mod_mpoly_vec_randtest_not_zero(fmpz_mod_mpoly_vec_t vec, flint_rand_t state, slong len, slong poly_len, ulong exp_bound, fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    fmpz_mod_mpoly_vec_set_length(vec, len, ctx);

    for (i = 0; i < len; i++)
    {
        do {
            fmpz_mod_mpoly_randtest_bound(vec->p + i, state, poly_len, exp_bound, ctx);
        } while (fmpz_mod_mpoly_is_zero(vec->p + i, ctx));
    }

    vec->length = len;
}

void
fmpz_mod_mpoly_vec_set_length(fmpz_mod_mpoly_vec_t vec, slong len, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    if (len > vec->length)
    {
        fmpz_mod_mpoly_vec_fit_length(vec, len, ctx);
        for (i = vec->length; i < len; i++)
            fmpz_mod_mpoly_zero(vec->p + i, ctx);
    }
    else if (len < vec->length)
    {
        for (i = len; i < vec->length; i++)
           fmpz_mod_mpoly_zero(vec->p + i, ctx);
    }

    vec->length = len;
}
