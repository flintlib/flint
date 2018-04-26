/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int _nmod_mpoly_divides_try_dense(slong * Adegs, slong * Bdegs, slong nvars,
                                                        slong Alen, slong Blen)
{
    slong i, total_dense_size;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    total_dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, total_dense_size, total_dense_size, Adegs[i] + 1);
        if (hi != WORD(0) || total_dense_size <= WORD(0))
            return 0;
    }

    return total_dense_size < WORD(5000000)
            && total_dense_size/Alen < WORD(10);
}

int nmod_mpoly_divides(nmod_mpoly_t Q,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int ret;
    slong i, * Adegs, * Bdegs;
    TMP_INIT;

    if (B->length == 0)
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in nmod_mpoly_divides");
    }

    TMP_START;

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS || A->length < 50)
    {
        ret = nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);

    } else
    {
        Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
        Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
        mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
        mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

        /* quick degree check */
        ret = -1;
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            if (Adegs[i] < Bdegs[i])
            {
                ret = 0;
                break;
            }
        }
        if (ret == -1)
        {
            if (_nmod_mpoly_divides_try_dense(Adegs, Bdegs, ctx->minfo->nvars,
                                                         A->length, B->length))
            {
                ret = nmod_mpoly_divides_dense(Q, A, B, ctx);
            }

            if (ret == -1)
            {
                ret = nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
            }
        }
    }

    TMP_END;
    return ret;
}

