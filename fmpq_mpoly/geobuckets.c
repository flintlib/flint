/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_geobucket_init(fmpq_mpoly_geobucket_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    B->length = 0;
}

void fmpq_mpoly_geobucket_clear(fmpq_mpoly_geobucket_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B->length; i++)
        fmpq_mpoly_clear(B->polys + i, ctx);
}

/* empty out bucket B into polynomial p */
void fmpq_mpoly_geobucket_empty(fmpq_mpoly_t p, fmpq_mpoly_geobucket_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpq_mpoly_zero(p, ctx);
    for (i = 0; i < B->length; i++)
    {
        fmpq_mpoly_add(p, p, B->polys + i, ctx);
        fmpq_mpoly_clear(B->polys + i, ctx);
    }
    B->length = 0;
}

void fmpq_mpoly_geobucket_print(fmpq_mpoly_geobucket_t B, const char ** x,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    printf("{");
    for (i = 0; i < B->length; i++)
    {
        fmpq_mpoly_print_pretty(B->polys + i, x, ctx);
        if (i + 1 < B->length)
            printf(", ");
    }
    printf("}");
}

/* ceiling(log_4(x)) - 1 */
slong fmpq_mpoly_geobucket_clog4(slong x)
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

void fmpq_mpoly_geobucket_fit_length(fmpq_mpoly_geobucket_t B, slong len,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong j;
    for (j = B->length; j < len; j++)
    {
        fmpq_mpoly_init(B->polys + j, ctx);
        fmpq_mpoly_zero(B->polys + j, ctx);
    }
    B->length = j;
}

/* set bucket B to polynomial p */
void fmpq_mpoly_geobucket_set(fmpq_mpoly_geobucket_t B, fmpq_mpoly_t p,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpq_mpoly_geobucket_clear(B, ctx);
    i = fmpq_mpoly_geobucket_clog4(p->zpoly->length);
    fmpq_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpq_mpoly_set(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
void _fmpq_mpoly_geobucket_fix(fmpq_mpoly_geobucket_t B, slong i,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    while (fmpq_mpoly_geobucket_clog4((B->polys + i)->zpoly->length) > i)
    {
        FLINT_ASSERT(i + 1 <= B->length);
        if (i + 1 == B->length)
        {
            fmpq_mpoly_init(B->polys + i + 1, ctx);
            fmpq_mpoly_zero(B->polys + i + 1, ctx);
            B->length = i + 2;
        }
        fmpq_mpoly_add(B->polys + i + 1, B->polys + i + 1, B->polys + i, ctx);
        fmpq_mpoly_zero(B->polys + i, ctx);
        i++;
    }
}

/* add polynomial p to buckect B */
void fmpq_mpoly_geobucket_add(fmpq_mpoly_geobucket_t B, fmpq_mpoly_t p,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    i = fmpq_mpoly_geobucket_clog4(p->zpoly->length);
    fmpq_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpq_mpoly_add(B->polys + i, B->polys + i, p, ctx);
    _fmpq_mpoly_geobucket_fix(B, i, ctx);
}

/* sub polynomial p to buckect B */
void fmpq_mpoly_geobucket_sub(fmpq_mpoly_geobucket_t B, fmpq_mpoly_t p,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    i = fmpq_mpoly_geobucket_clog4(p->zpoly->length);
    fmpq_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fmpq_mpoly_sub(B->polys + i, B->polys + i, p, ctx);
    _fmpq_mpoly_geobucket_fix(B, i, ctx);
}

void fmpq_mpoly_geobucket_set_fmpz(fmpq_mpoly_geobucket_t B, fmpz_t c,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fmpq_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fmpq_mpoly_clear(B->polys + i, ctx);
    B->length = 1;
    fmpq_mpoly_set_fmpz(B->polys + 0, c, ctx);

}

void fmpq_mpoly_geobucket_gen(fmpq_mpoly_geobucket_t B, slong var,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fmpq_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fmpq_mpoly_clear(B->polys + i, ctx);
    B->length = 1;
    fmpq_mpoly_gen(B->polys + 0, var, ctx);
}

void fmpq_mpoly_geobucket_add_inplace(fmpq_mpoly_geobucket_t B1,
                         fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B2->length; i++)
        fmpq_mpoly_geobucket_add(B1, B2->polys + i, ctx);

}

void fmpq_mpoly_geobucket_sub_inplace(fmpq_mpoly_geobucket_t B1,
                         fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B2->length; i++)
        fmpq_mpoly_geobucket_sub(B1, B2->polys + i, ctx);
}

void fmpq_mpoly_geobucket_neg_inplace(fmpq_mpoly_geobucket_t B1,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B1->length; i++)
        fmpq_mpoly_neg(B1->polys + i, B1->polys + i, ctx);
}

void fmpq_mpoly_geobucket_mul_inplace(fmpq_mpoly_geobucket_t B1,
                         fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_mpoly_t a, b;
    fmpq_mpoly_init(a, ctx);
    fmpq_mpoly_init(b, ctx);

    fmpq_mpoly_geobucket_empty(a, B1, ctx);
    fmpq_mpoly_geobucket_empty(b, B2, ctx);
    fmpq_mpoly_mul(a, a, b, ctx);
    fmpq_mpoly_geobucket_set(B1, a, ctx);

    fmpq_mpoly_clear(a, ctx);
    fmpq_mpoly_clear(b, ctx);
}

void fmpq_mpoly_geobucket_pow_fmpz_inplace(fmpq_mpoly_geobucket_t B1,
                                    const fmpz_t k, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_mpoly_t a;
    fmpq_mpoly_init(a, ctx);

    fmpq_mpoly_geobucket_empty(a, B1, ctx);
    fmpq_mpoly_pow_fmpz(a, a, k, ctx);
    fmpq_mpoly_geobucket_set(B1, a, ctx);

    fmpq_mpoly_clear(a, ctx);
}


int fmpq_mpoly_geobucket_divides_inplace(fmpq_mpoly_geobucket_t B1,
                         fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx)
{
    int ret;
    fmpq_mpoly_t a, b;
    fmpq_mpoly_init(a, ctx);
    fmpq_mpoly_init(b, ctx);

    fmpq_mpoly_geobucket_empty(a, B1, ctx);
    fmpq_mpoly_geobucket_empty(b, B2, ctx);

    ret = fmpq_mpoly_divides(a, a, b, ctx);
    fmpq_mpoly_geobucket_set(B1, a, ctx);

    fmpq_mpoly_clear(a, ctx);
    fmpq_mpoly_clear(b, ctx);
    return ret;
}
