/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


/*
    unpack the field-wise minimum of poly_exps into max_fields
    must have len > 0
*/
void mpoly_min_fields_ui(ulong * min_fields, const ulong * poly_exps,
                           slong len, mp_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmin, mask;
    TMP_INIT;

    FLINT_ASSERT(len > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp(bits, mctx);

    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
    {
        mask = (mask << bits) + (UWORD(1) << (bits - 1));
    }

    TMP_START;

    pmin = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    i = 0;
    mpoly_monomial_set(pmin, poly_exps + N*i, N);
    for (i = 1; i < len; i++)
    {
        mpoly_monomial_min(pmin, pmin, poly_exps + N*i, bits, N, mask);
    }

    mpoly_unpack_vec_ui(min_fields, pmin, bits, mctx->nfields, 1);

    TMP_END;
}


void mpoly_min_fields_fmpz(fmpz * min_fields, const ulong * poly_exps,
                           slong len, mp_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, j, N;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    if (bits <= FLINT_BITS)
    {
        ulong * min_uis = (ulong *) TMP_ALLOC(mctx->nfields*sizeof(ulong));
        mpoly_min_fields_ui(min_uis, poly_exps, len, bits, mctx);
        for (j = 0; j < mctx->nfields; j++)
        {
            fmpz_set_ui(min_fields + j, min_uis[j]);
        }
    }
    else
    {
        fmpz * tmp_exps;
        tmp_exps = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
        for (j = 0; j < mctx->nfields; j++)
        {
            fmpz_init(tmp_exps + j);
        }

        N = mpoly_words_per_exp(bits, mctx);

        i = 0;
        mpoly_unpack_vec_fmpz(min_fields, poly_exps + N*i, bits, mctx->nfields, 1);
        for (i = 1; i < len; i++)
        {
            mpoly_unpack_vec_fmpz(tmp_exps, poly_exps + N*i, bits, mctx->nfields, 1);
            for (j = 0; j < mctx->nfields; j++)
            {
                if (fmpz_cmp(min_fields + j, tmp_exps + j) > 0)
                    fmpz_set(min_fields + j, tmp_exps + j);
            }
        }

        for (j = 0; j < mctx->nfields; j++)
        {
            fmpz_clear(tmp_exps + j);
        }
    }

    TMP_END;
}
