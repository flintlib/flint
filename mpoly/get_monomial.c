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
#include "mpoly.h"


void mpoly_get_monomial_ui(ulong * user_exps, const ulong * poly_exps,
                                            slong bits, const mpoly_ctx_t mctx)
{
    slong i;
    ulong * tmp_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(mctx->nfields*sizeof(ulong));
    mpoly_unpack_vec_ui(tmp_exps, poly_exps, bits, mctx->nfields, 1);

    for (i = 0; i < mctx->nvars; i++)
        user_exps[i] = tmp_exps[mctx->rev ? i : mctx->nvars - 1 - i];

    TMP_END;
}


void mpoly_get_monomial_fmpz(fmpz * user_exps, const ulong * poly_exps,
                                      mp_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, nvars = mctx->nvars;
    fmpz * tmp_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(tmp_exps + i);

    mpoly_unpack_vec_fmpz(tmp_exps, poly_exps, bits, mctx->nfields, 1);

    for (i = 0; i < nvars; i++)
        fmpz_set(user_exps + i, tmp_exps + (mctx->rev ? i : nvars - 1 - i));

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(tmp_exps + i);
    TMP_END;
}
