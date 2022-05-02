/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

/* unpack the field-wise maximum of poly_exps into max_fields */
void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps,
                        slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmax, mask;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, mctx);

    mask = mpoly_overflow_mask_sp(bits);

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
                        slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmax, mask;
    TMP_INIT;

    TMP_START;

    N = mpoly_words_per_exp(bits, mctx);
    pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (i = 0; i < N; i++)
        pmax[i] = 0;

    if (bits <= FLINT_BITS)
    {
        mask = mpoly_overflow_mask_sp(bits);

        for (i = 0; i < len; i++)
            mpoly_monomial_max(pmax, pmax, poly_exps + N*i, bits, N, mask);
    }
    else
    {
        for (i = 0; i < len; i++)
            mpoly_monomial_max_mp(pmax, pmax, poly_exps + N*i, bits, N);
    }

    mpoly_unpack_vec_fmpz(max_fields, pmax, bits, mctx->nfields, 1);

    TMP_END;
}
