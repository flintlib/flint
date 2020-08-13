/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

int mpoly_is_gen(ulong * exps, slong var, flint_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    int ret;
    slong i;
    fmpz * unpacked_exps;
    TMP_INIT;

    TMP_START;

    unpacked_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(unpacked_exps + i);

    mpoly_get_monomial_ffmpz(unpacked_exps, exps, bits, mctx);

    if (var >= 0)
    {
        ret = 1;
        for (i = 0; i < mctx->nvars; i++)
        {
            if (!fmpz_equal_si(unpacked_exps + i, i == var))
            {
                ret = 0;
                break;
            }
        }
    }
    else
    {
        int count = 0;
        for (i = 0; i < mctx->nvars; i++)
        {
            if (fmpz_is_one(unpacked_exps + i))
            {
                count++;
                if (count > 1)
                    break;
            }
            else if (!fmpz_is_zero(unpacked_exps + i))
            {
                count = 2;
                break;
            }
        }
        ret = (count == 1);
    }

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(unpacked_exps + i);

    TMP_END;
    return ret;
}
