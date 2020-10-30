/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


slong mpoly_degree_si(const ulong * exps,
                slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx)
{
    if (len == 0)
        return -WORD(1);

    /* sometimes we don't have to look very far */
    if (mctx->ord == ORD_LEX && var == 0)
        len = 1;

    if (bits <= FLINT_BITS)
    {
        slong r;
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        slong offset, shift, N, i;

        mpoly_gen_offset_shift_sp(&offset, &shift, var, bits, mctx);
        N = mpoly_words_per_exp_sp(bits, mctx);

        i = 0;
        r = (exps[N*i + offset] >> shift) & mask;
        for (i = 1; i < len; i++)
        {
            ulong k = (exps[N*i + offset] >> shift) & mask;
            if (r < k)
                r = k;
        }
        return r;
    }
    else
    {
        slong * degs, r;
        TMP_INIT;
        TMP_START;
        degs = (slong *) TMP_ALLOC(mctx->nvars*sizeof(slong));
        mpoly_degrees_si(degs, exps, len, bits, mctx);
        r = degs[var];
        TMP_END;
        return r;
    }
}


void mpoly_degree_fmpz(fmpz_t deg, const ulong * exps,
                slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * degs;
    TMP_INIT;

    TMP_START;
    degs = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(degs + i);

    mpoly_degrees_ffmpz(degs, exps, len, bits, mctx);

    fmpz_swap(deg, degs + var);
    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(degs + i);

    TMP_END;
}
