/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

/* get the offset and shift where the variable of index var is stored */
void mpoly_gen_offset_shift(slong * offset, slong * shift,
                  slong idx, slong N, mp_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong fpw = FLINT_BITS/bits;

    FLINT_ASSERT(bits <= FLINT_BITS);

    if (!mctx->rev)
        idx = nvars - 1 - idx;

    *offset = idx/fpw;
    *shift  = idx%fpw*bits;
}

/* additionally get the monomial as well */
void mpoly_gen_oneexp_offset_shift(ulong * oneexp, slong * offset, slong * shift,
                  slong idx, slong N, mp_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong fpw = FLINT_BITS/bits;
    slong i;

    FLINT_ASSERT(bits <= FLINT_BITS);

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

/* just get the monomial */
void mpoly_gen_monomial_sp(ulong * oneexp, slong var, mp_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    ulong nvars = mctx->nvars;
    ulong fpw = FLINT_BITS/bits;
    ulong idx = var;
    slong i, N;

    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, mctx);
    for (i = 0; i < N; i++)
        oneexp[i] = 0;

    if (!mctx->rev)
        idx = nvars - 1 - idx;

    oneexp[idx/fpw] |= UWORD(1) << (idx%fpw*bits);
    if (mctx->deg)
        oneexp[nvars/fpw] |= UWORD(1) << (nvars%fpw*bits);
}

/* get the offset where the variable of index var is stored */
slong mpoly_gen_offset_mp(slong idx, slong N, mp_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong wpf = bits/FLINT_BITS;

    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(bits % FLINT_BITS == WORD(0));

    if (!mctx->rev)
        idx = nvars - 1 - idx;
    return idx*wpf;
}

/* additionally get the monomial as well */
void mpoly_gen_oneexp_offset_mp(ulong * oneexp, slong * offset,
                  slong idx, slong N, mp_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong wpf = bits/FLINT_BITS;
    slong i;

    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(bits % FLINT_BITS == WORD(0));

    for (i = 0; i < N; i++)
        oneexp[i] = 0;

    if (!mctx->rev)
        idx = nvars - 1 - idx;

    *offset = idx*wpf;

    oneexp[idx*wpf] = UWORD(1);
    if (mctx->deg)
        oneexp[nvars*wpf] = UWORD(1);
}

/* just get the monomial */
void mpoly_gen_monomial_mp(ulong * oneexp, slong var, mp_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    ulong nvars = mctx->nvars;
    ulong wpf = bits/FLINT_BITS;
    ulong idx;
    slong i, N;

    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(bits % FLINT_BITS == WORD(0));

    N = mpoly_words_per_exp_mp(bits, mctx);
    for (i = 0; i < N; i++)
        oneexp[i] = 0;

    idx = var;
    if (!mctx->rev)
        idx = nvars - 1 - var;

    oneexp[idx*wpf] = UWORD(1);
    if (mctx->deg)
        oneexp[nvars*wpf] = UWORD(1);
}

