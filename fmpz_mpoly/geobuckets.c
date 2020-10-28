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

/* ceiling(log_4(x)) - 1 */
slong fmpz_mpoly_geobucket_clog4(slong x)
{
    if (x <= 4)
        return 0;
    /*
        FLINT_BIT_COUNT returns unsigned int.
        Signed division is not defined.
        Do the calculation with unsigned ints and then convert to slong.
    */
    x = (FLINT_BIT_COUNT(x - 1) - (unsigned int) 1)/((unsigned int) 2);
    return x;
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
    i = fmpz_mpoly_geobucket_clog4(p->length);
    B->length = 0;
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_swap(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
static void _fmpz_mpoly_geobucket_fix(fmpz_mpoly_geobucket_t B, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    while (fmpz_mpoly_geobucket_clog4((B->polys + i)->length) > i)
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

    i = fmpz_mpoly_geobucket_clog4(p->length);
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

    i = fmpz_mpoly_geobucket_clog4(p->length);
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_sub(B->temps + i, B->polys + i, p, ctx);
    fmpz_mpoly_swap(B->polys + i, B->temps + i, ctx);
    _fmpz_mpoly_geobucket_fix(B, i, ctx);
}

void fmpz_mpoly_geobucket_set_fmpz(fmpz_mpoly_geobucket_t B, fmpz_t c,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    B->length = 1;
    fmpz_mpoly_set_fmpz(B->polys + 0, c, ctx);
}

void fmpz_mpoly_geobucket_gen(fmpz_mpoly_geobucket_t B, slong var,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    B->length = 1;
    fmpz_mpoly_gen(B->polys + 0, var, ctx);
}

void fmpz_mpoly_geobucket_add_inplace(fmpz_mpoly_geobucket_t B1,
                         fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B2->length; i++)
        fmpz_mpoly_geobucket_add(B1, B2->polys + i, ctx);
}

void fmpz_mpoly_geobucket_sub_inplace(fmpz_mpoly_geobucket_t B1,
                         fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B2->length; i++)
        fmpz_mpoly_geobucket_sub(B1, B2->polys + i, ctx);
}

void fmpz_mpoly_geobucket_neg_inplace(fmpz_mpoly_geobucket_t B1,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B1->length; i++)
        fmpz_mpoly_neg(B1->polys + i, B1->polys + i, ctx);
}

void fmpz_mpoly_geobucket_mul_inplace(fmpz_mpoly_geobucket_t B1,
                         fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t a, b;
    fmpz_mpoly_init(a, ctx);
    fmpz_mpoly_init(b, ctx);

    fmpz_mpoly_geobucket_empty(a, B1, ctx);
    fmpz_mpoly_geobucket_empty(b, B2, ctx);
    fmpz_mpoly_mul(a, a, b, ctx);
    fmpz_mpoly_geobucket_set(B1, a, ctx);

    fmpz_mpoly_clear(a, ctx);
    fmpz_mpoly_clear(b, ctx);
}

void fmpz_mpoly_geobucket_pow_ui_inplace(fmpz_mpoly_geobucket_t B1,
                                           ulong k, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t a;
    fmpz_mpoly_init(a, ctx);

    fmpz_mpoly_geobucket_empty(a, B1, ctx);
    if (!fmpz_mpoly_pow_ui(a, a, k, ctx))
        flint_throw(FLINT_ERROR, "fmpz_mpoly_pow_ui failed");
    fmpz_mpoly_geobucket_set(B1, a, ctx);

    fmpz_mpoly_clear(a, ctx);
}

void fmpz_mpoly_geobucket_pow_fmpz_inplace(fmpz_mpoly_geobucket_t B1,
                                    const fmpz_t k, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t a;
    fmpz_mpoly_init(a, ctx);

    fmpz_mpoly_geobucket_empty(a, B1, ctx);
    if (!fmpz_mpoly_pow_fmpz(a, a, k, ctx))
        flint_throw(FLINT_ERROR, "fmpz_mpoly_pow_fmpz failed");
    fmpz_mpoly_geobucket_set(B1, a, ctx);

    fmpz_mpoly_clear(a, ctx);
}


int fmpz_mpoly_geobucket_divides_inplace(fmpz_mpoly_geobucket_t B1,
                         fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx)
{
    int ret;
    fmpz_mpoly_t a, b;
    fmpz_mpoly_init(a, ctx);
    fmpz_mpoly_init(b, ctx);

    fmpz_mpoly_geobucket_empty(a, B1, ctx);
    fmpz_mpoly_geobucket_empty(b, B2, ctx);

    ret = fmpz_mpoly_divides_monagan_pearce(a, a, b, ctx);
    fmpz_mpoly_geobucket_set(B1, a, ctx);

    fmpz_mpoly_clear(a, ctx);
    fmpz_mpoly_clear(b, ctx);
    return ret;
}

