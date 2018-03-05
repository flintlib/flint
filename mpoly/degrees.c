/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"


void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmax, mask;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            user_degs[i] = -WORD(1);
        return;
    }

    TMP_START;

    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    N = mpoly_words_per_exp(bits, mctx);
    pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (i = 0; i < N; i++)
        pmax[i] = 0;
    for (i = 0; i < len; i++)
        mpoly_monomial_max(pmax, pmax, poly_exps + N*i, bits, N, mask);

    mpoly_get_monomial_ui((ulong *) user_degs, pmax, bits, mctx);

    TMP_END;
}


void mpoly_degrees_fmpz(fmpz * user_degs, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i, j, N;
    fmpz * tmp_exps;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            fmpz_set_si(user_degs + i, -WORD(1));
        return;
    }

    TMP_START;
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (j = 0; j < mctx->nvars; j++) {
        fmpz_zero(user_degs + j);
        fmpz_init(tmp_exps + j);
    }

    N = mpoly_words_per_exp(bits, mctx);

    for (i = 0; i < len; i++)
    {
        mpoly_get_monomial_fmpz(tmp_exps, poly_exps + N*i, bits, mctx);
        for (j = 0; j < mctx->nvars; j++)
            if (fmpz_cmp(user_degs + j, tmp_exps + j) < 0)
                fmpz_set(user_degs + j, tmp_exps + j);
    }

    for (j = 0; j < mctx->nvars; j++)
        fmpz_clear(tmp_exps + j);


    TMP_END;
}
