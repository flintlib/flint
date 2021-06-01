/*
    Copyright (C) 2018-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* TODO: this could a bit faster */
int mpoly_is_poly(
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    slong var,
    const mpoly_ctx_t mctx)
{
    int ret = 1;
    slong i, j;
    slong N = mpoly_words_per_exp(Abits, mctx);
    slong nvars = mctx->nvars;
    fmpz * t;
    TMP_INIT;

    TMP_START;
    t = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init(t + i);

    for (i = 0; i < Alen; i++)
    {
        mpoly_get_monomial_ffmpz(t, Aexps + N*i, Abits, mctx);
        for (j = 0; j < nvars; j++)
        {
            if (j != var && !fmpz_is_zero(t + j))
            {
                ret = 0;
                goto cleanup;
            }
        }
    }

cleanup:

    for (i = 0; i < nvars; i++)
        fmpz_clear(t + i);

    TMP_END;

    return ret;
}
