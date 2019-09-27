/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpoly_geobucket_init(fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    B->length = 0;
}

void fq_nmod_mpoly_geobucket_clear(fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B->length; i++)
        fq_nmod_mpoly_clear(B->polys + i, ctx);
}

/* empty out bucket B into polynomial p */
void fq_nmod_mpoly_geobucket_empty(fq_nmod_mpoly_t p, fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpoly_zero(p, ctx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_add(p, p, B->polys + i, ctx);
        fq_nmod_mpoly_clear(B->polys + i, ctx);
    }
    B->length = 0;
}

void fq_nmod_mpoly_geobucket_print(fq_nmod_mpoly_geobucket_t B, const char ** x,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    printf("{");
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_print_pretty(B->polys + i, x, ctx);
        if (i + 1 < B->length)
            printf(", ");
    }
    printf("}");
}

/* ceiling(log_4(x)) - 1 */
slong fq_nmod_mpoly_geobucket_clog4(slong x)
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

void fq_nmod_mpoly_geobucket_fit_length(fq_nmod_mpoly_geobucket_t B, slong len,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong j;
    for (j = B->length; j < len; j++)
    {
        fq_nmod_mpoly_init(B->polys + j, ctx);
        fq_nmod_mpoly_zero(B->polys + j, ctx);
    }
    B->length = j;
}

/* set bucket B to polynomial p */
void fq_nmod_mpoly_geobucket_set(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpoly_geobucket_clear(B, ctx);
    i = fq_nmod_mpoly_geobucket_clog4(p->length);
    fq_nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fq_nmod_mpoly_set(B->polys + i, p, ctx);
    B->length = i + 1;
}

/* internal function for fixing overflows */
void _fq_nmod_mpoly_geobucket_fix(fq_nmod_mpoly_geobucket_t B, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    while (fq_nmod_mpoly_geobucket_clog4((B->polys + i)->length) > i)
    {
        FLINT_ASSERT(i + 1 <= B->length);
        if (i + 1 == B->length)
        {
            fq_nmod_mpoly_init(B->polys + i + 1, ctx);
            fq_nmod_mpoly_zero(B->polys + i + 1, ctx);
            B->length = i + 2;
        }
        fq_nmod_mpoly_add(B->polys + i + 1, B->polys + i + 1, B->polys + i, ctx);
        fq_nmod_mpoly_zero(B->polys + i, ctx);
        i++;
    }
}

/* add polynomial p to buckect B */
void fq_nmod_mpoly_geobucket_add(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    i = fq_nmod_mpoly_geobucket_clog4(p->length);
    fq_nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fq_nmod_mpoly_add(B->polys + i, B->polys + i, p, ctx);
    _fq_nmod_mpoly_geobucket_fix(B, i, ctx);
}

/* sub polynomial p to buckect B */
void fq_nmod_mpoly_geobucket_sub(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    i = fq_nmod_mpoly_geobucket_clog4(p->length);
    fq_nmod_mpoly_geobucket_fit_length(B, i + 1, ctx);
    fq_nmod_mpoly_sub(B->polys + i, B->polys + i, p, ctx);
    _fq_nmod_mpoly_geobucket_fix(B, i, ctx);
}

void fq_nmod_mpoly_geobucket_set_ui(fq_nmod_mpoly_geobucket_t B, ulong c,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fq_nmod_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fq_nmod_mpoly_clear(B->polys + i, ctx);
    B->length = 1;
    fq_nmod_mpoly_set_ui(B->polys + 0, c, ctx);
}

void fq_nmod_mpoly_geobucket_set_fq_nmod_gen(fq_nmod_mpoly_geobucket_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fq_nmod_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fq_nmod_mpoly_clear(B->polys + i, ctx);
    B->length = 1;
    fq_nmod_mpoly_set_fq_nmod_gen(B->polys + 0, ctx);
}

void fq_nmod_mpoly_geobucket_gen(fq_nmod_mpoly_geobucket_t B, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (B->length == 0)
        fq_nmod_mpoly_init(B->polys + 0, ctx);
    for (i = 1; i < B->length; i++)
        fq_nmod_mpoly_clear(B->polys + i, ctx);
    B->length = 1;
    fq_nmod_mpoly_gen(B->polys + 0, var, ctx);
}

void fq_nmod_mpoly_geobucket_add_inplace(fq_nmod_mpoly_geobucket_t A,
                    fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B->length; i++)
        fq_nmod_mpoly_geobucket_add(A, B->polys + i, ctx);

}

void fq_nmod_mpoly_geobucket_sub_inplace(fq_nmod_mpoly_geobucket_t A,
                    fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < B->length; i++)
        fq_nmod_mpoly_geobucket_sub(A, B->polys + i, ctx);
}

void fq_nmod_mpoly_geobucket_neg_inplace(fq_nmod_mpoly_geobucket_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fq_nmod_mpoly_neg(A->polys + i, A->polys + i, ctx);
}

void fq_nmod_mpoly_geobucket_mul_inplace(fq_nmod_mpoly_geobucket_t A,
                    fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_t a, b;
    fq_nmod_mpoly_init(a, ctx);
    fq_nmod_mpoly_init(b, ctx);

    fq_nmod_mpoly_geobucket_empty(a, A, ctx);
    fq_nmod_mpoly_geobucket_empty(b, B, ctx);
    fq_nmod_mpoly_mul(a, a, b, ctx);
    fq_nmod_mpoly_geobucket_set(A, a, ctx);

    fq_nmod_mpoly_clear(a, ctx);
    fq_nmod_mpoly_clear(b, ctx);
}

void fq_nmod_mpoly_geobucket_pow_ui(fq_nmod_mpoly_geobucket_t A,
                                        ulong k, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_t a;
    fq_nmod_mpoly_init(a, ctx);

    fq_nmod_mpoly_geobucket_empty(a, A, ctx);
    fq_nmod_mpoly_pow_ui(a, a, k, ctx);
    fq_nmod_mpoly_geobucket_set(A, a, ctx);

    fq_nmod_mpoly_clear(a, ctx);
}

void fq_nmod_mpoly_geobucket_pow_fmpz_inplace(fq_nmod_mpoly_geobucket_t A,
                                 const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_t a;
    fq_nmod_mpoly_init(a, ctx);

    fq_nmod_mpoly_geobucket_empty(a, A, ctx);
    fq_nmod_mpoly_pow_fmpz(a, a, k, ctx);
    fq_nmod_mpoly_geobucket_set(A, a, ctx);

    fq_nmod_mpoly_clear(a, ctx);
}


int fq_nmod_mpoly_geobucket_divides_inplace(fq_nmod_mpoly_geobucket_t A,
                    fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    int ret = 0;
    fq_nmod_mpoly_t a, b;
    fq_nmod_mpoly_init(a, ctx);
    fq_nmod_mpoly_init(b, ctx);

    fq_nmod_mpoly_geobucket_empty(a, A, ctx);
    fq_nmod_mpoly_geobucket_empty(b, B, ctx);

    if (!fq_nmod_mpoly_is_zero(b, ctx))
    {
        ret = fq_nmod_mpoly_divides(a, a, b, ctx);
        fq_nmod_mpoly_geobucket_set(A, a, ctx);
    }

    fq_nmod_mpoly_clear(a, ctx);
    fq_nmod_mpoly_clear(b, ctx);
    return ret;
}
