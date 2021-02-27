/*
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpoly_geobucket_init(fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < FLINT_BITS/2; i++)
    {
        fq_nmod_mpoly_init(B->polys + i, ctx);
        fq_nmod_mpoly_init(B->temps + i, ctx);
    }
    B->length = 0;
}

void fq_nmod_mpoly_geobucket_clear(fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < FLINT_BITS/2; i++)
    {
        fq_nmod_mpoly_clear(B->polys + i, ctx);
        fq_nmod_mpoly_clear(B->temps + i, ctx);
    }
}

/* empty out bucket B into polynomial p */
void fq_nmod_mpoly_geobucket_empty(fq_nmod_mpoly_t p, fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (B->length < 2)
    {
        if (B->length < 1)
            fq_nmod_mpoly_zero(p, ctx);
        else
            fq_nmod_mpoly_set(p, B->polys + 0, ctx);
    }
    else if (B->length == 2)
    {
        fq_nmod_mpoly_add(p, B->polys + 1, B->polys + 0, ctx);
    }
    else
    {
        fq_nmod_mpoly_add(B->temps + 1, B->polys + 1, B->polys + 0, ctx);
        for (i = 2; i < B->length - 1; i++)
            fq_nmod_mpoly_add(B->temps + i, B->polys + i, B->temps + i - 1, ctx);
        fq_nmod_mpoly_add(p, B->polys + i, B->temps + i - 1, ctx);
    }

    B->length = 0;
}

void fq_nmod_mpoly_geobucket_fit_length(fq_nmod_mpoly_geobucket_t B, slong len,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong j;
    for (j = B->length; j < len; j++)
        fq_nmod_mpoly_zero(B->polys + j, ctx);
    B->length = j;
}

/* set bucket B to polynomial p */
void fq_nmod_mpoly_geobucket_set(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i = mpoly_geobucket_clog4(p->length);
    B->length = 0;
    fq_nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fq_nmod_mpoly_swap(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
void _fq_nmod_mpoly_geobucket_fix(fq_nmod_mpoly_geobucket_t B, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    while (mpoly_geobucket_clog4((B->polys + i)->length) > i)
    {
        FLINT_ASSERT(i + 1 <= B->length);
        if (i + 1 == B->length)
        {
            B->length = i + 2;
            fq_nmod_mpoly_set(B->polys + i + 1, B->polys + i, ctx);
        }
        else
        {
            fq_nmod_mpoly_add(B->temps + i + 1, B->polys + i + 1, B->polys + i, ctx);
            fq_nmod_mpoly_swap(B->polys + i + 1, B->temps + i + 1, ctx);
        }
        fq_nmod_mpoly_zero(B->polys + i, ctx);
        i++;
    }
}

/* add polynomial p to buckect B */
void fq_nmod_mpoly_geobucket_add(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (fq_nmod_mpoly_is_zero(p, ctx))
        return;

    i = mpoly_geobucket_clog4(p->length);
    fq_nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fq_nmod_mpoly_add(B->temps + i, B->polys + i, p, ctx);
    fq_nmod_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _fq_nmod_mpoly_geobucket_fix(B, i, ctx);
}

/* sub polynomial p to buckect B */
void fq_nmod_mpoly_geobucket_sub(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (fq_nmod_mpoly_is_zero(p, ctx))
        return;

    i = mpoly_geobucket_clog4(p->length);
    fq_nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fq_nmod_mpoly_sub(B->temps + i, B->polys + i, p, ctx);
    fq_nmod_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _fq_nmod_mpoly_geobucket_fix(B, i, ctx);
}

