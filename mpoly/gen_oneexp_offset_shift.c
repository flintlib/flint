/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_gen_offset_shift(slong * offset, slong * shift,
                       slong idx, slong N, slong bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong fpw = FLINT_BITS/bits;

    if (!mctx->rev)
        idx = nvars - 1 - idx;

    *offset = idx/fpw;
    *shift  = idx%fpw*bits;
}



void mpoly_gen_oneexp_offset_shift(ulong * oneexp, slong * offset, slong * shift,
                       slong idx, slong N, slong bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong fpw = FLINT_BITS/bits;
    slong i;

    for (i = 0; i < N; i++)
        oneexp[i] = 0;

    if (!mctx->rev)
        idx = nvars - 1 - idx;

    *offset = idx/fpw;
    *shift  = idx%fpw*bits;

    oneexp[idx/fpw] |= UWORD(1) << (idx%fpw*bits);
    if (mctx->deg)
        oneexp[nvars/fpw] |= UWORD(1) << (nvars%fpw*bits);
}
