/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_gen_offset_shift(slong * _offset, slong * _shift,
                       slong idx, slong N, slong bits, const mpoly_ctx_t mctx)
{
    slong nfields = mctx->nfields;
    int deg = mctx->deg;
    int rev = mctx->rev;
    slong fpw = FLINT_BITS/bits;
    slong offset, shift;

    if (rev)
    {
        offset = (nfields - 1 - idx)/fpw;
        shift  = (nfields - 1 - idx)%fpw;
    } else
    {
        offset = (deg + idx)/fpw;
        shift  = (deg + idx)%fpw;
    }
    shift = (fpw - 1 - shift) * bits;

    * _offset = offset;
    * _shift = shift;
}



void mpoly_gen_oneexp_offset_shift(ulong * oneexp, slong * _offset, slong * _shift,
                       slong idx, slong N, slong bits, const mpoly_ctx_t mctx)
{
    slong nfields = mctx->nfields;
    int deg = mctx->deg;
    int rev = mctx->rev;
    slong fpw = FLINT_BITS/bits;
    slong i, offset, shift;

    if (rev)
    {
        offset = (nfields - 1 - idx)/fpw;
        shift  = (nfields - 1 - idx)%fpw;
    } else
    {
        offset = (deg + idx)/fpw;
        shift  = (deg + idx)%fpw;
    }
    shift = (fpw - 1 - shift) * bits;

    * _offset = offset;
    * _shift = shift;

    for (i = 0; i < N; i++)
        oneexp[i] = 0;
    oneexp[offset] = WORD(1) << shift;
    if (deg)
        oneexp[0] |= WORD(1) << ((fpw - 1)*bits);
}
