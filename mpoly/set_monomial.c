/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

void mpoly_set_monomial_ui(ulong * poly_exps, const ulong * user_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong nfields = mctx->nfields;
    slong i = 0;
    ulong * tmp_exps, degree;
    fmpz * big_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

    degree = 0;
    for (i = 0; i < nvars; i++) {
        degree += user_exps[i];
        if (mctx->deg && degree < user_exps[i])
            goto big_case;
        tmp_exps[mctx->rev ? i : nvars - 1 - i] = user_exps[i];
    }

    if (mctx->deg)
        tmp_exps[nvars] = degree;

    mpoly_pack_vec_ui(poly_exps, tmp_exps, bits, nfields, 1);

done:

    TMP_END;
    return;


big_case:

    big_exps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));

    for (i = 0; i < nvars; i++)
        fmpz_init_set_ui(big_exps + i, user_exps[i]);

    mpoly_set_monomial_ffmpz(poly_exps, big_exps, bits, mctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(big_exps + i);

    goto done;
}


void mpoly_set_monomial_ffmpz(ulong * poly_exps, const fmpz * user_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong nfields = mctx->nfields;
    slong i = 0;
    fmpz * tmp_exps;
    fmpz_t degree;
    TMP_INIT;

    TMP_START;
    fmpz_init_set_ui(degree, 0);
    tmp_exps = (fmpz *) TMP_ALLOC(nfields*sizeof(fmpz));
    for (i = 0; i < nvars; i++) {
        fmpz_add(degree, degree, user_exps + i);
        fmpz_init_set(tmp_exps + (mctx->rev ? i : nvars - 1 - i), user_exps + i);
    }

    if (mctx->deg)
        fmpz_init_set(tmp_exps + nvars, degree);

    mpoly_pack_vec_fmpz(poly_exps, tmp_exps, bits, nfields, 1);

    fmpz_clear(degree);
    for (i = 0; i < nvars; i++)
        fmpz_clear(tmp_exps + i);
    if (mctx->deg)
        fmpz_clear(tmp_exps + nvars);

    TMP_END;
}

void mpoly_set_monomial_pfmpz(ulong * poly_exps, fmpz * const * user_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong nfields = mctx->nfields;
    slong i = 0;
    fmpz * tmp_exps;
    fmpz_t degree;
    TMP_INIT;

    TMP_START;
    fmpz_init_set_ui(degree, 0);
    tmp_exps = (fmpz *) TMP_ALLOC(nfields*sizeof(fmpz));
    for (i = 0; i < nvars; i++) {
        fmpz_add(degree, degree, user_exps[i]);
        fmpz_init_set(tmp_exps + (mctx->rev ? i : nvars - 1 - i), user_exps[i]);
    }

    if (mctx->deg)
        fmpz_init_set(tmp_exps + nvars, degree);

    mpoly_pack_vec_fmpz(poly_exps, tmp_exps, bits, nfields, 1);

    fmpz_clear(degree);
    for (i = 0; i < nvars; i++)
        fmpz_clear(tmp_exps + i);
    if (mctx->deg)
        fmpz_clear(tmp_exps + nvars);

    TMP_END;
}
