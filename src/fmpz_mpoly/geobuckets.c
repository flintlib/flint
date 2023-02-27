/*
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"


void fmpz_mpoly_geobucket_init(fmpz_mpoly_geobucket_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < FLINT_BITS/2; i++)
    {
        fmpz_mpoly_init(B->polys + i, ctx);
        fmpz_mpoly_init(B->temps + i, ctx);
    }
    B->length = 0;
}

void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < FLINT_BITS/2; i++)
    {
        fmpz_mpoly_clear(B->polys + i, ctx);
        fmpz_mpoly_clear(B->temps + i, ctx);
    }
}

/* empty out bucket B into polynomial p */
void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p, fmpz_mpoly_geobucket_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (B->length < 2)
    {
        if (B->length < 1)
            fmpz_mpoly_zero(p, ctx);
        else
            fmpz_mpoly_set(p, B->polys + 0, ctx);
    }
    else if (B->length == 2)
    {
        fmpz_mpoly_add(p, B->polys + 1, B->polys + 0, ctx);
    }
    else
    {
        fmpz_mpoly_add(B->temps + 1, B->polys + 1, B->polys + 0, ctx);
        for (i = 2; i < B->length - 1; i++)
            fmpz_mpoly_add(B->temps + i, B->polys + i, B->temps + i - 1, ctx);
        fmpz_mpoly_add(p, B->polys + i, B->temps + i - 1, ctx);
    }

    B->length = 0;
}

void fmpz_mpoly_geobucket_fit_length(fmpz_mpoly_geobucket_t B, slong len,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong j;
    for (j = B->length; j < len; j++)
        fmpz_mpoly_zero(B->polys + j, ctx);
    B->length = j;
}

/* set bucket B to polynomial p */
void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    i = mpoly_geobucket_clog4(p->length);
    B->length = 0;
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_swap(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
static void _fmpz_mpoly_geobucket_fix(fmpz_mpoly_geobucket_t B, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    while (mpoly_geobucket_clog4((B->polys + i)->length) > i)
    {
        FLINT_ASSERT(i + 1 <= B->length);
        if (i + 1 == B->length)
        {
            B->length = i + 2;
            fmpz_mpoly_set(B->polys + i + 1, B->polys + i, ctx);
        }
        else
        {
            fmpz_mpoly_add(B->temps + i + 1, B->polys + i + 1, B->polys + i, ctx);
            fmpz_mpoly_swap(B->polys + i + 1, B->temps + i + 1, ctx);
        }
        fmpz_mpoly_zero(B->polys + i, ctx);
        i++;
    }
}

/* add polynomial p to buckect B */
void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (p->length < 1)
        return;

    i = mpoly_geobucket_clog4(p->length);
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_add(B->temps + i, B->polys + i, p, ctx);
    fmpz_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _fmpz_mpoly_geobucket_fix(B, i, ctx);
}

/* sub polynomial p to buckect B */
void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (p->length < 1)
        return;

    i = mpoly_geobucket_clog4(p->length);
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_sub(B->temps + i, B->polys + i, p, ctx);
    fmpz_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _fmpz_mpoly_geobucket_fix(B, i, ctx);
}

