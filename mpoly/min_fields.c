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

/* unpack the field-wise minimum of poly_exps into min_fields */
void mpoly_min_fields_ui_sp(ulong * min_fields, const ulong * poly_exps,
                           slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmin, mask;
    TMP_INIT;

    FLINT_ASSERT(len > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, mctx);

    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    TMP_START;

    pmin = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_set(pmin, poly_exps + N*0, N);
    for (i = 1; i < len; i++)
        mpoly_monomial_min(pmin, pmin, poly_exps + N*i, bits, N, mask);

    mpoly_unpack_vec_ui(min_fields, pmin, bits, mctx->nfields, 1);

    TMP_END;
}


void mpoly_min_fields_fmpz(fmpz * min_fields, const ulong * poly_exps,
                           slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    ulong * pmin, mask;
    TMP_INIT;

    FLINT_ASSERT(len > 0);

    TMP_START;

    N = mpoly_words_per_exp(bits, mctx);
    pmin = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_set(pmin, poly_exps + N*0, N);

    if (bits <= FLINT_BITS)
    {
        mask = 0;
        for (i = 0; i < FLINT_BITS/bits; i++)
            mask = (mask << bits) + (UWORD(1) << (bits - 1));

        for (i = 1; i < len; i++)
            mpoly_monomial_min(pmin, pmin, poly_exps + N*i, bits, N, mask);
    }
    else
    {
        for (i = 1; i < len; i++)
            mpoly_monomial_min_mp(pmin, pmin, poly_exps + N*i, bits, N);
    }

    mpoly_unpack_vec_fmpz(min_fields, pmin, bits, mctx->nfields, 1);

    TMP_END;
}
