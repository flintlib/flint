/*
    Copyright (C) 2018 Daniel Schultz
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


slong mpoly_totaldegree_si(const ulong * exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong r;
    fmpz_t td;
    fmpz_init(td);
    mpoly_totaldegree_fmpz(td, exps, len, bits, mctx);
    r = fmpz_get_si(td);
    fmpz_clear(td);
    return r;
}

void mpoly_totaldegree_fmpz(fmpz_t totdeg, const ulong * exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i, j, N;
    fmpz_t tot;
    fmpz * tmp_exps;
    TMP_INIT;
    fmpz_set_si(totdeg, -WORD(1));
    TMP_START;
    fmpz_init(tot);
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (j = 0; j < mctx->nvars; j++)
        fmpz_init(tmp_exps + j);
    N = mpoly_words_per_exp(bits, mctx);
    for (i = 0; i < len; i++)
    {
        mpoly_get_monomial_ffmpz(tmp_exps, exps + N*i, bits, mctx);
        fmpz_zero(tot);
        for (j = 0; j < mctx->nvars; j++)
            fmpz_add(tot, tot, tmp_exps + j);
        if (fmpz_cmp(totdeg, tot) < 0)
            fmpz_swap(totdeg, tot);
    }
    fmpz_clear(tot);
    for (j = 0; j < mctx->nvars; j++)
        fmpz_clear(tmp_exps + j);
    TMP_END;
}
