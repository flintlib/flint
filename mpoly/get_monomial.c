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


void mpoly_get_monomial(ulong * user_exps, const ulong * poly_exps,
                                            slong bits, const mpoly_ctx_t mctx)
{
    slong i;
    ulong * tmp_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(mctx->nfields*sizeof(ulong));
    mpoly_unpack_vec(tmp_exps, poly_exps, bits, mctx->nfields, 1);

    if (mctx->rev)
    {
        for (i = mctx->nfields - 1; i >= mctx->deg; i--)
            user_exps[mctx->nfields - i - 1] = tmp_exps[i];
    } else
    {
        for (i = mctx->deg; i < mctx->nfields; i++)
            user_exps[i - mctx->deg] = tmp_exps[i];
    }

    TMP_END;
}
