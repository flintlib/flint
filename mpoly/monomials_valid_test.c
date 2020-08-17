/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


int mpoly_monomials_valid_test(ulong * exps, slong len,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    int ret =  1;
    slong N, i, j;
    fmpz * fields;
    TMP_INIT;

    if (!mctx->deg)
        return ret;

    TMP_START;
    fields = TMP_ALLOC(mctx->nfields*sizeof(fmpz));

    for (j = 0; j < mctx->nfields; j++)
        fmpz_init(fields + j);

    N = mpoly_words_per_exp(bits, mctx);
    for (i = 0; i < len; i++)
    {
        mpoly_unpack_vec_fmpz(fields, exps + i*N, bits, mctx->nfields, 1);

        for (j = 0; j < mctx->nvars; j++)
        {
            fmpz_sub(fields + mctx->nvars, fields + mctx->nvars, fields + j);
        }

        if (!fmpz_is_zero(fields + mctx->nvars))
        {
            ret = 0;
            break;
        }
    }

    for (j = 0; j < mctx->nfields; j++)
        fmpz_clear(fields + j);

    TMP_END;
    return ret;
}
