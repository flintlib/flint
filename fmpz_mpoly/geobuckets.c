/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
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
    B->length = 0;
}

void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B->length; i++)
        fmpz_mpoly_clear(B->polys + i, ctx);
}

/* empty out bucket B into polynomial p */
void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p, fmpz_mpoly_geobucket_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_zero(p, ctx);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_add(p, p, B->polys + i, ctx);
        fmpz_mpoly_clear(B->polys + i, ctx);
    }
    B->length = 0;
}

void fmpz_mpoly_geobucket_print(fmpz_mpoly_geobucket_t B, const char ** x,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    printf("{");
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_print_pretty(B->polys + i, x, ctx);
        if (i + 1 < B->length)
            printf(", ");
    }
    printf("}");
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
    {
        fmpz_mpoly_init(B->polys + j, ctx);
        fmpz_mpoly_zero(B->polys + j, ctx);
    }
    B->length = j;
}

/* set bucket B to polynomial p */
void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_geobucket_clear(B, ctx);
    i = fmpz_mpoly_geobucket_clog4(p->length);
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_set(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
void _fmpz_mpoly_geobucket_fix(fmpz_mpoly_geobucket_t B, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    while (fmpz_mpoly_geobucket_clog4((B->polys + i)->length) > i)
    {
        FLINT_ASSERT(i + 1 <= B->length);
        if (i + 1 == B->length)
        {
            fmpz_mpoly_init(B->polys + i + 1, ctx);
            fmpz_mpoly_zero(B->polys + i + 1, ctx);
            B->length = i + 2;
        }
        fmpz_mpoly_add(B->polys + i + 1, B->polys + i + 1, B->polys + i, ctx);
        fmpz_mpoly_zero(B->polys + i, ctx);
        i++;
    }
}

/* add polynomial p to buckect B */
void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    i = fmpz_mpoly_geobucket_clog4(p->length);
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_add(B->polys + i, B->polys + i, p, ctx);
    _fmpz_mpoly_geobucket_fix(B, i, ctx);
}

/* sub polynomial p to buckect B */
void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    i = fmpz_mpoly_geobucket_clog4(p->length);
    fmpz_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpz_mpoly_sub(B->polys + i, B->polys + i, p, ctx);
    _fmpz_mpoly_geobucket_fix(B, i, ctx);
}

void fmpz_mpoly_geobucket_set_fmpz(fmpz_mpoly_geobucket_t B, fmpz_t c,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fmpz_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fmpz_mpoly_clear(B->polys + i, ctx);
    B->length = 1;
    fmpz_mpoly_set_fmpz(B->polys + 0, c, ctx);
}

void fmpz_mpoly_geobucket_gen(fmpz_mpoly_geobucket_t B, slong var,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fmpz_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fmpz_mpoly_clear(B->polys + i, ctx);
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
    fmpz_mpoly_mul_johnson(a, a, b, ctx);
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
    fmpz_mpoly_pow_fps(a, a, k, ctx);
    fmpz_mpoly_geobucket_set(B1, a, ctx);

    fmpz_mpoly_clear(a, ctx);
}

void fmpz_mpoly_geobucket_pow_fmpz_inplace(fmpz_mpoly_geobucket_t B1,
                                    const fmpz_t k, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t a;
    fmpz_mpoly_init(a, ctx);

    fmpz_mpoly_geobucket_empty(a, B1, ctx);
    fmpz_mpoly_pow_fmpz(a, a, k, ctx);
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
