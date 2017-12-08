/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


void mpoly_get_cmpmask(ulong * cmpmask, slong N, slong bits,
                                                        const mpoly_ctx_t mctx)
{
    slong i;

    /* LEX and DEGLEX use a clear mask */
    i = 0;
    do {
        cmpmask[i] = WORD(0);
    } while (++i < N);

    /* DEGREVLEX has a bit set everywhere except the in most significant field */
    if (mctx->ord == ORD_DEGREVLEX)
    {
        if (bits <= FLINT_BITS)
        {
            cmpmask[0] = (UWORD(1) << ((bits)*((FLINT_BITS)/(bits) - 1))) - 1;
            for (i = 1; i < N; i++)
                cmpmask[i] = -UWORD(1);
        } else {

            flint_throw(FLINT_ERROR, "bits > FLINT_BITS in mpoly_get_cmpmask");

            for (i = 0; i < N - bits/FLINT_BITS; i++)
                cmpmask[i] = -UWORD(1);
        }
    }
}
