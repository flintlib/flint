/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef UTILS_FLINT_H
#define UTILS_FLINT_H

#ifdef UTILS_FLINT_INLINES_C
#define UTILS_FLINT_INLINE
#else
#define UTILS_FLINT_INLINE static __inline__
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "flint/fmpz_mpoly.h"
#include "flint/fmpq.h"
#include "calcium.h"

/* Multivariate polynomials */

void fmpz_mpoly_symmetric_gens(fmpz_mpoly_t res, ulong k, slong * vars, slong n, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_symmetric(fmpz_mpoly_t res, ulong k, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_spoly(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx);

/* Vectors of multivariate polynomials */

typedef struct
{
    fmpz_mpoly_struct * p;
    slong length;
    slong alloc;
}
fmpz_mpoly_vec_struct;

typedef fmpz_mpoly_vec_struct fmpz_mpoly_vec_t[1];

#define fmpz_mpoly_vec_entry(vec, i) ((vec)->p + (i))

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_init(fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)
{
    vec->p = NULL;
    vec->length = 0;
    vec->alloc = 0;
}

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_print(const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    flint_printf("[");
    for (i = 0; i < F->length; i++)
    {
        fmpz_mpoly_print_pretty(F->p + i, NULL, ctx);
        if (i < F->length - 1)
            flint_printf(", ");
    }
    flint_printf("]");
}

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_swap(fmpz_mpoly_vec_t x, fmpz_mpoly_vec_t y, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_t tmp;
    *tmp = *x;
    *x = *y;
    *y = *tmp;
}

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_fit_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)
{
    if (len > vec->alloc)
    {
        slong i;

        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->p = flint_realloc(vec->p, len * sizeof(fmpz_mpoly_struct));

        for (i = vec->alloc; i < len; i++)
            fmpz_mpoly_init(vec->p + i, ctx);

        vec->alloc = len;
    }
}

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_clear(fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < vec->alloc; i++)
        fmpz_mpoly_clear(vec->p + i, ctx);

    flint_free(vec->p);
}

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_set(fmpz_mpoly_vec_t dest, const fmpz_mpoly_vec_t src, const fmpz_mpoly_ctx_t ctx)
{
    if (dest != src)
    {
        slong i;

        fmpz_mpoly_vec_fit_length(dest, src->length, ctx);

        for (i = 0; i < src->length; i++)
            fmpz_mpoly_set(dest->p + i, src->p + i, ctx);

        dest->length = src->length;
    }
}

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_append(fmpz_mpoly_vec_t vec, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_fit_length(vec, vec->length + 1, ctx);
    fmpz_mpoly_set(vec->p + vec->length, f, ctx);
    vec->length++;
}

UTILS_FLINT_INLINE slong
fmpz_mpoly_vec_insert_unique(fmpz_mpoly_vec_t vec, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < vec->length; i++)
    {
        if (fmpz_mpoly_equal(vec->p + i, f, ctx))
            return i;
    }

    fmpz_mpoly_vec_append(vec, f, ctx);
    return vec->length - 1;
}

void fmpz_mpoly_vec_set_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx);

UTILS_FLINT_INLINE void
fmpz_mpoly_vec_randtest_not_zero(fmpz_mpoly_vec_t vec, flint_rand_t state, slong len, slong poly_len, slong bits, ulong exp_bound, fmpz_mpoly_ctx_t ctx)
{
    slong i;

    fmpz_mpoly_vec_set_length(vec, len, ctx);

    for (i = 0; i < len; i++)
    {
        do {
            fmpz_mpoly_randtest_bound(vec->p + i, state, poly_len, bits, exp_bound, ctx);
        } while (fmpz_mpoly_is_zero(vec->p + i, ctx));
    }

    vec->length = len;
}

void fmpz_mpoly_vec_set_primitive_unique(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

/* Index pairs (for Buchberger algorithm) */

typedef struct
{
    slong a;
    slong b;
} pair_t;

typedef struct
{
    pair_t * pairs;
    slong length;
    slong alloc;
}
pairs_struct;

typedef pairs_struct pairs_t[1];

UTILS_FLINT_INLINE void
pairs_init(pairs_t vec)
{
    vec->pairs = NULL;
    vec->length = 0;
    vec->alloc = 0;
}

UTILS_FLINT_INLINE void
pairs_fit_length(pairs_t vec, slong len)
{
    if (len > vec->alloc)
    {
        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->pairs = flint_realloc(vec->pairs, len * sizeof(pair_t));
        vec->alloc = len;
    }
}

UTILS_FLINT_INLINE void
pairs_clear(pairs_t vec)
{
    flint_free(vec->pairs);
}

UTILS_FLINT_INLINE void
pairs_append(pairs_t vec, slong i, slong j)
{
    pairs_fit_length(vec, vec->length + 1);
    vec->pairs[vec->length].a = i;
    vec->pairs[vec->length].b = j;
    vec->length++;
}

UTILS_FLINT_INLINE void
pairs_insert_unique(pairs_t vec, slong i, slong j)
{
    slong k;

    for (k = 0; k < vec->length; k++)
    {
        if (vec->pairs[k].a == i && vec->pairs[k].b == j)
            return;
    }

    pairs_append(vec, i, j);
}

/* Ideals and Groebner bases */

void fmpz_mpoly_reduction_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_vec_t I, const fmpz_mpoly_ctx_t ctx);
int fmpz_mpoly_vec_is_groebner(const fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);
pair_t fmpz_mpoly_select_pop_pair(pairs_t pairs, const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);
int fmpz_mpoly_buchberger_naive_with_limits(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F,
    slong ideal_len_limit, slong poly_len_limit, slong poly_bits_limit, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_vec_autoreduction(fmpz_mpoly_vec_t H, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_vec_autoreduction_groebner(fmpz_mpoly_vec_t H, const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx);

int fmpz_mpoly_vec_is_autoreduced(const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

