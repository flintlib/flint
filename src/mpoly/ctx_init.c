/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_ctx_init(mpoly_ctx_t mctx, slong nvars, const ordering_t ord)
{
    flint_bitcnt_t bits;

    mctx->ord = ord;
    if (nvars > 0)
    {
        mctx->nvars = nvars;
        switch (ord)
        {
            case ORD_LEX:
                mctx->deg = 0;
                mctx->rev = 0;
                break;
            case ORD_DEGLEX:
                mctx->deg = 1;
                mctx->rev = 0;
                break;
            case ORD_DEGREVLEX:
                mctx->deg = 1;
                mctx->rev = 1;
                break;
            default:
                flint_throw(FLINT_ERROR, "Invalid ordering in mpoly_ctx_init");
        }
    }
    else
    {
        mctx->nvars = 0;
        mctx->deg = 1;  /* degree field will always be zero */
        mctx->rev = 0;
    }
    mctx->nfields = mctx->nvars + mctx->deg;

    for (bits = 1; bits <= FLINT_BITS; bits++)
    {
        mctx->lut_words_per_exp[bits - 1]
                                   = (mctx->nfields - 1)/(FLINT_BITS/bits) + 1;
    }

    for (bits = 1; bits <= FLINT_BITS; bits++)
    {
        flint_bitcnt_t new_bits = FLINT_MAX(bits, MPOLY_MIN_BITS);
        while (new_bits < FLINT_BITS && mctx->lut_words_per_exp[new_bits - 1]
                                       == mctx->lut_words_per_exp[new_bits])
        {
            new_bits++;
        }
        mctx->lut_fix_bits[bits - 1] = new_bits;
    }

}
