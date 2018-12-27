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

/* unpack the field-wise maximum of poly_exps into max_fields */
void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmax, mask;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, mctx);

    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    TMP_START;

    pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (i = 0; i < N; i++)
        pmax[i] = 0;
    for (i = 0; i < len; i++)
        mpoly_monomial_max(pmax, pmax, poly_exps + N*i, bits, N, mask);

    mpoly_unpack_vec_ui(max_fields, pmax, bits, mctx->nfields, 1);

    TMP_END;
}


void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps,
                                 slong len, slong bits, const mpoly_ctx_t mctx)
{
    slong i, j, N;
    ulong * pmax, mask;
    fmpz * tmp_exps;
    TMP_INIT;

    TMP_START;

    if (bits <= FLINT_BITS)
    {
        N = mpoly_words_per_exp_sp(bits, mctx);

        mask = 0;
        for (i = 0; i < FLINT_BITS/bits; i++)
            mask = (mask << bits) + (UWORD(1) << (bits - 1));

        pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
        for (i = 0; i < N; i++)
            pmax[i] = 0;
        for (i = 0; i < len; i++)
            mpoly_monomial_max(pmax, pmax, poly_exps + N*i, bits, N, mask);

        mpoly_unpack_vec_fmpz(max_fields, pmax, bits, mctx->nfields, 1);
    }
    else
    {
        N = mpoly_words_per_exp_mp(bits, mctx);

        tmp_exps = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
        for (j = 0; j < mctx->nfields; j++)
        {
            fmpz_zero(max_fields + j);
            fmpz_init(tmp_exps + j);
        }

        for (i = 0; i < len; i++)
        {
            mpoly_unpack_vec_fmpz(tmp_exps, poly_exps + N*i, bits, mctx->nfields, 1);
            for (j = 0; j < mctx->nfields; j++)
            {
                if (fmpz_cmp(max_fields + j, tmp_exps + j) < 0)
                    fmpz_set(max_fields + j, tmp_exps + j);
            }
        }

        for (j = 0; j < mctx->nfields; j++)
        {
            fmpz_clear(tmp_exps + j);
        }
    }

    TMP_END;
}
