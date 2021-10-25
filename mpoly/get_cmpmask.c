/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_get_cmpmask(ulong * cmpmask, slong N, flint_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    slong i;

    if (mctx->rev)
    {
        /* DEGREVLEX has a bit set everywhere except the in most significant field */
        if (bits <= FLINT_BITS)
        {
            for (i = 0; i + 1 < N; i++)
                cmpmask[i] = -UWORD(1);
            cmpmask[N - 1] = (UWORD(1) << (mctx->nvars%(FLINT_BITS/bits)*bits)) - UWORD(1);
        }
        else
        {
            for (i = 0; i < N - bits/FLINT_BITS; i++)
                cmpmask[i] = -UWORD(1);
            for (; i < N; i++)
                cmpmask[i] = UWORD(0);
        }
    }
    else
    {
        /* LEX and DEGLEX use a clear mask */
        for (i = 0; i < N; i++)
            cmpmask[i] = UWORD(0);
    }
}
