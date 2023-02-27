/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

int mpoly_degrees_fit_si(const ulong * poly_exps, slong len,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, j, N;
    int ret;
    fmpz * tmp_exps;
    TMP_INIT;

    if (len == 0)
    {
        return 1;
    }

    TMP_START;
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (j = 0; j < mctx->nvars; j++) {
        fmpz_init(tmp_exps + j);
    }

    N = mpoly_words_per_exp(bits, mctx);

    ret = 1;
    for (i = 0; i < len; i++)
    {
        mpoly_get_monomial_ffmpz(tmp_exps, poly_exps + N*i, bits, mctx);
        for (j = 0; j < mctx->nvars; j++)
        {
            if (!fmpz_fits_si(tmp_exps + j))
            {
                ret = 0;
                break;
            }
        }
    }

    for (j = 0; j < mctx->nvars; j++)
        fmpz_clear(tmp_exps + j);

    TMP_END;
    return ret;
}
