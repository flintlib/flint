/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

/* get the offset and shift where variable var is stored in packed form */
void mpoly_gen_offset_shift_sp(slong * offset, slong * shift, slong var,
                                      flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    ulong fpw = FLINT_BITS/bits;
    ulong idx;

    FLINT_ASSERT(0 <= var && var < nvars);
    FLINT_ASSERT(bits <= FLINT_BITS);

    idx = var;
    if (!mctx->rev)
        idx = nvars - 1 - var;

    *offset = idx/fpw;
    *shift  = idx%fpw*bits;
}

/* additionally get the monomial as well */
void mpoly_gen_monomial_offset_shift_sp(ulong * mexp, slong * offset,
            slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    ulong idx;
    ulong nvars = mctx->nvars;
    ulong fpw = FLINT_BITS/bits;
    slong i, N;

    FLINT_ASSERT(0 <= var && var < nvars);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, mctx);
    for (i = 0; i < N; i++)
        mexp[i] = 0;

    idx = var;
    if (!mctx->rev)
        idx = nvars - 1 - var;

    *offset = idx/fpw;
    *shift  = idx%fpw*bits;

    mexp[idx/fpw] |= UWORD(1) << (idx%fpw*bits);
    if (mctx->deg)
        mexp[nvars/fpw] |= UWORD(1) << (nvars%fpw*bits);
}

/* just get the monomial */
void mpoly_gen_monomial_sp(ulong * mexp, slong var, flint_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    ulong nvars = mctx->nvars;
    ulong fpw = FLINT_BITS/bits;
    ulong idx = var;
    slong i, N;

    FLINT_ASSERT(0 <= var && var < nvars);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, mctx);
    for (i = 0; i < N; i++)
        mexp[i] = 0;

    idx = var;
    if (!mctx->rev)
        idx = nvars - 1 - var;

    mexp[idx/fpw] |= UWORD(1) << (idx%fpw*bits);
    if (mctx->deg)
        mexp[nvars/fpw] |= UWORD(1) << (nvars%fpw*bits);
}

/* get the offset where the variable of index var is stored */
slong mpoly_gen_offset_mp(slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    ulong wpf = bits/FLINT_BITS;
    ulong idx;

    FLINT_ASSERT(0 <= var && var < nvars);
    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(bits % FLINT_BITS == WORD(0));

    idx = var;
    if (!mctx->rev)
        idx = nvars - 1 - idx;

    return idx*wpf;
}

/* additionally get the monomial as well */
slong mpoly_gen_monomial_offset_mp(ulong * mexp,
                           slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    ulong nvars = mctx->nvars;
    ulong wpf = bits/FLINT_BITS;
    ulong idx;
    slong i, N, offset;

    FLINT_ASSERT(0 <= var && var < nvars);
    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(bits % FLINT_BITS == WORD(0));

    N = mpoly_words_per_exp_mp(bits, mctx);
    for (i = 0; i < N; i++)
        mexp[i] = 0;

    idx = var;
    if (!mctx->rev)
        idx = nvars - 1 - var;

    offset = idx*wpf;

    mexp[idx*wpf] = UWORD(1);
    if (mctx->deg)
        mexp[nvars*wpf] = UWORD(1);

    return offset;
}
