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


void mpoly_degrees(slong * user_degs, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmax, mask;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nfields - mctx->deg; i++)
            user_degs[i] = -WORD(1);
        return;
    }

    TMP_START;

    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    N = words_per_exp(mctx->nfields, bits);
    pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (i = 0; i < N; i++)
        pmax[i] = 0;
    for (i = 0; i < len; i++)
        mpoly_monomial_max(pmax, pmax, poly_exps + N*i, bits, N, mask);

    mpoly_get_monomial((ulong *) user_degs, pmax, bits, mctx);

    TMP_END;
}
