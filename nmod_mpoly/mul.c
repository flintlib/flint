/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int _nmod_mpoly_mul_try_dense(slong * Adegs, slong * Bdegs,
                            slong Alen, slong Blen, const nmod_mpoly_ctx_t ctx)
{
    slong i, total_dense_size;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    total_dense_size = WORD(1);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        umul_ppmm(hi, total_dense_size, total_dense_size, Adegs[i] + Bdegs[i] + 1);
        if (hi != WORD(0) || total_dense_size <= WORD(0))
            return 0;
    }

    return total_dense_size < WORD(5000000)
            && total_dense_size/Alen/Blen < WORD(10);
}


int _nmod_mpoly_mul_try_array(slong * Adegs, slong * Bdegs,
                            slong Alen, slong Blen, const nmod_mpoly_ctx_t ctx)
{
    slong i, total_dense_size, total_degree;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    if (ctx->minfo->nvars <= WORD(1) || ctx->minfo->nvars >= WORD(8))
    {
        return 0;
    }

    /* accept array method if the array is probably at least 10% full */

    if (ctx->minfo->ord == ORD_LEX)
    {
        total_dense_size = WORD(1);
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            umul_ppmm(hi, total_dense_size, total_dense_size, Adegs[i] + Bdegs[i] + 1);
            if (hi != WORD(0) || total_dense_size <= WORD(0))
                return 0;
        }
        return total_dense_size <= WORD(5000000)
                && total_dense_size/Alen/Blen < WORD(10);

    } else if (ctx->minfo->ord == ORD_DEGLEX || ctx->minfo->ord == ORD_DEGREVLEX)
    {
        total_degree = 0;
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            total_degree += Adegs[i];
            if (total_degree < 0)
                return 0;
            total_degree += Bdegs[i];
            if (total_degree < 0)
                return 0;
        }

        /* the relevent portion of the array has approx size d^nvars/nvars!*/
        total_dense_size = WORD(1);
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            umul_ppmm(hi, total_dense_size, total_dense_size, total_degree);
            if (hi != WORD(0) || total_dense_size < WORD(0))
                return 0;
        }
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            total_dense_size /= i+1;
        }

        return total_dense_size <= WORD(5000000)
                && total_dense_size/Alen/Blen < WORD(10);

    } else
    {
        return 0;
    }
}

void nmod_mpoly_mul(nmod_mpoly_t P,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * Adegs, * Bdegs;
    TMP_INIT;

    TMP_START;

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS
             || A->length < 100 || B->length < 100)
    {
        nmod_mpoly_mul_johnson(P, A, B, ctx);

    } else
    {
        Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
        Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
        mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
        mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

        success = 0;
        if (_nmod_mpoly_mul_try_dense(Adegs, Bdegs, A->length, B->length, ctx))
        {
            success = nmod_mpoly_mul_dense(P, A, B, ctx);
        }

        if (!success)
        {
            if (_nmod_mpoly_mul_try_array(Adegs, Bdegs, A->length, B->length, ctx))
            {
                success = nmod_mpoly_mul_array(P, A, B, ctx);
            }

            if (!success)
            {
                nmod_mpoly_mul_johnson(P, A, B, ctx);
            }
        }
    }

    TMP_END;
}

