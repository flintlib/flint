/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


void mpoly_set_monomial_ui(ulong * poly_exps, const ulong * user_exps,
                                            slong bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong nfields = mctx->nfields;
    slong i = 0;
    ulong * tmp_exps, degree;
    TMP_INIT;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

    degree = 0;
    for (i = 0; i < nvars; i++) {
        degree += user_exps[i];
        tmp_exps[mctx->rev ? i : nvars - 1 - i] = user_exps[i];
    }

    if (mctx->deg)
        tmp_exps[nvars] = degree;

    mpoly_pack_vec_ui(poly_exps, tmp_exps, bits, nfields, 1);

    TMP_END;
}
