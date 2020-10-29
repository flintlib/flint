/*
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void nmod_mpoly_geobucket_init(nmod_mpoly_geobucket_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < FLINT_BITS/2; i++)
    {
        nmod_mpoly_init(B->polys + i, ctx);
        nmod_mpoly_init(B->temps + i, ctx);
    }
    B->length = 0;
}

void nmod_mpoly_geobucket_clear(nmod_mpoly_geobucket_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < FLINT_BITS/2; i++)
    {
        nmod_mpoly_clear(B->polys + i, ctx);
        nmod_mpoly_clear(B->temps + i, ctx);
    }
}

/* empty out bucket B into polynomial p */
void nmod_mpoly_geobucket_empty(nmod_mpoly_t p, nmod_mpoly_geobucket_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (B->length < 2)
    {
        if (B->length < 1)
            nmod_mpoly_zero(p, ctx);
        else
            nmod_mpoly_set(p, B->polys + 0, ctx);
    }
    else if (B->length == 2)
    {
        nmod_mpoly_add(p, B->polys + 1, B->polys + 0, ctx);
    }
    else
    {
        nmod_mpoly_add(B->temps + 1, B->polys + 1, B->polys + 0, ctx);
        for (i = 2; i < B->length - 1; i++)
            nmod_mpoly_add(B->temps + i, B->polys + i, B->temps + i - 1, ctx);
        nmod_mpoly_add(p, B->polys + i, B->temps + i - 1, ctx);
    }

    B->length = 0;
}

void nmod_mpoly_geobucket_fit_length(nmod_mpoly_geobucket_t B, slong len,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong j;
    for (j = B->length; j < len; j++)
        nmod_mpoly_zero(B->polys + j, ctx);
    B->length = j;
}

/* set bucket B to polynomial p */
void nmod_mpoly_geobucket_set(nmod_mpoly_geobucket_t B, nmod_mpoly_t p,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i = mpoly_geobucket_clog4(p->length);
    B->length = 0;
    nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    nmod_mpoly_swap(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
void _nmod_mpoly_geobucket_fix(nmod_mpoly_geobucket_t B, slong i,
                                                    const nmod_mpoly_ctx_t ctx)
{
    while (mpoly_geobucket_clog4((B->polys + i)->length) > i)
    {
        FLINT_ASSERT(i + 1 <= B->length);
        if (i + 1 == B->length)
        {
            B->length = i + 2;
            nmod_mpoly_set(B->polys + i + 1, B->polys + i, ctx);
        }
        else
        {
            nmod_mpoly_add(B->temps + i + 1, B->polys + i + 1, B->polys + i, ctx);
            nmod_mpoly_swap(B->polys + i + 1, B->temps + i + 1, ctx);
        }
        nmod_mpoly_zero(B->polys + i, ctx);
        i++;
    }
}

/* add polynomial p to bucket B */
void nmod_mpoly_geobucket_add(nmod_mpoly_geobucket_t B, nmod_mpoly_t p,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (nmod_mpoly_is_zero(p, ctx))
        return;

    i = mpoly_geobucket_clog4(p->length);
    nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    nmod_mpoly_add(B->temps + i, B->polys + i, p, ctx);
    nmod_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _nmod_mpoly_geobucket_fix(B, i, ctx);
}

/* sub polynomial p to buckect B */
void nmod_mpoly_geobucket_sub(nmod_mpoly_geobucket_t B, nmod_mpoly_t p,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (nmod_mpoly_is_zero(p, ctx))
        return;

    i = mpoly_geobucket_clog4(p->length);
    nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    nmod_mpoly_sub(B->temps + i, B->polys + i, p, ctx);
    nmod_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _nmod_mpoly_geobucket_fix(B, i, ctx);
}

void nmod_mpoly_geobucket_set_ui(nmod_mpoly_geobucket_t B, ulong c,
                                                    const nmod_mpoly_ctx_t ctx)
{
    B->length = 1;
    nmod_mpoly_set_ui(B->polys + 0, c, ctx);
}

void nmod_mpoly_geobucket_gen(nmod_mpoly_geobucket_t B, slong var,
                                                    const nmod_mpoly_ctx_t ctx)
{
    B->length = 1;
    nmod_mpoly_gen(B->polys + 0, var, ctx);
}

void nmod_mpoly_geobucket_add_inplace(nmod_mpoly_geobucket_t B1,
                         nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B2->length; i++)
        nmod_mpoly_geobucket_add(B1, B2->polys + i, ctx);
}

void nmod_mpoly_geobucket_sub_inplace(nmod_mpoly_geobucket_t B1,
                         nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B2->length; i++)
        nmod_mpoly_geobucket_sub(B1, B2->polys + i, ctx);
}

void nmod_mpoly_geobucket_neg_inplace(nmod_mpoly_geobucket_t B1,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B1->length; i++)
        nmod_mpoly_neg(B1->polys + i, B1->polys + i, ctx);
}

void nmod_mpoly_geobucket_mul_inplace(nmod_mpoly_geobucket_t B1,
                         nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t a, b;
    nmod_mpoly_init(a, ctx);
    nmod_mpoly_init(b, ctx);

    nmod_mpoly_geobucket_empty(a, B1, ctx);
    nmod_mpoly_geobucket_empty(b, B2, ctx);
    nmod_mpoly_mul_johnson(a, a, b, ctx);
    nmod_mpoly_geobucket_set(B1, a, ctx);

    nmod_mpoly_clear(a, ctx);
    nmod_mpoly_clear(b, ctx);
}

void nmod_mpoly_geobucket_pow_ui(nmod_mpoly_geobucket_t B1,
                                           ulong k, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t a;
    nmod_mpoly_init(a, ctx);

    nmod_mpoly_geobucket_empty(a, B1, ctx);
    if (!nmod_mpoly_pow_ui(a, a, k, ctx))
        flint_throw(FLINT_ERROR, "nmod_mpoly_pow_ui failed");
    nmod_mpoly_geobucket_set(B1, a, ctx);

    nmod_mpoly_clear(a, ctx);
}

void nmod_mpoly_geobucket_pow_fmpz_inplace(nmod_mpoly_geobucket_t B1,
                                    const fmpz_t k, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t a;
    nmod_mpoly_init(a, ctx);

    nmod_mpoly_geobucket_empty(a, B1, ctx);
    if (!nmod_mpoly_pow_fmpz(a, a, k, ctx))
        flint_throw(FLINT_ERROR, "nmod_mpoly_pow_fmpz failed");
    nmod_mpoly_geobucket_set(B1, a, ctx);

    nmod_mpoly_clear(a, ctx);
}


int nmod_mpoly_geobucket_divides_inplace(nmod_mpoly_geobucket_t B1,
                         nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx)
{
    int ret = 0;
    nmod_mpoly_t a, b;
    nmod_mpoly_init(a, ctx);
    nmod_mpoly_init(b, ctx);

    nmod_mpoly_geobucket_empty(a, B1, ctx);
    nmod_mpoly_geobucket_empty(b, B2, ctx);

    nmod_mpoly_mul_johnson(a, a, b, ctx);
    nmod_mpoly_geobucket_set(B1, a, ctx);

    nmod_mpoly_clear(a, ctx);
    nmod_mpoly_clear(b, ctx);
    return ret;
}
