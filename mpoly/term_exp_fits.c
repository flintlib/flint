/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

int mpoly_term_exp_fits_ui(ulong * exps, flint_bitcnt_t bits,
                                               slong n, const mpoly_ctx_t mctx)
{
    slong i, N;
    int ret;
    fmpz * unpacked_exps;
    TMP_INIT;

    TMP_START;

    unpacked_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(unpacked_exps + i);

    N = mpoly_words_per_exp(bits, mctx);
    mpoly_get_monomial_ffmpz(unpacked_exps, exps + N*n, bits, mctx);

    ret = 1;
    for (i = 0; i < mctx->nvars; i++)
    {
        ret = ret && fmpz_abs_fits_ui(unpacked_exps + i);
        fmpz_clear(unpacked_exps + i);
    }

    TMP_END;
    return ret;
}

int mpoly_term_exp_fits_si(ulong * exps, flint_bitcnt_t bits,
                                               slong n, const mpoly_ctx_t mctx)
{
    slong i, N;
    int ret;
    fmpz * unpacked_exps;
    TMP_INIT;

    TMP_START;

    unpacked_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(unpacked_exps + i);

    N = mpoly_words_per_exp(bits, mctx);
    mpoly_get_monomial_ffmpz(unpacked_exps, exps + N*n, bits, mctx);

    ret = 1;
    for (i = 0; i < mctx->nvars; i++)
    {
        ret = ret && fmpz_fits_si(unpacked_exps + i);
        fmpz_clear(unpacked_exps + i);
    }

    TMP_END;
    return ret;
}
