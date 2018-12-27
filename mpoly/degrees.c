/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * max_fields;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            user_degs[i] = -WORD(1);
        return;
    }

    TMP_START;

    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields + i);

    mpoly_max_fields_fmpz(max_fields, poly_exps, len, bits, mctx);

    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)user_degs, max_fields, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields + i);

    TMP_END;
}


void mpoly_degrees_ffmpz(fmpz * user_degs, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * max_fields;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            fmpz_set_si(user_degs + i, -WORD(1));
        return;
    }

    TMP_START;

    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields + i);

    mpoly_max_fields_fmpz(max_fields, poly_exps, len, bits, mctx);

    mpoly_get_monomial_ffmpz_unpacked_ffmpz(user_degs, max_fields, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields + i);

    TMP_END;
}

void mpoly_degrees_pfmpz(fmpz ** user_degs, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * max_fields;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            fmpz_set_si(user_degs[i], -WORD(1));
        return;
    }

    TMP_START;

    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields + i);

    mpoly_max_fields_fmpz(max_fields, poly_exps, len, bits, mctx);

    mpoly_get_monomial_pfmpz_unpacked_ffmpz(user_degs, max_fields, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields + i);

    TMP_END;
}
