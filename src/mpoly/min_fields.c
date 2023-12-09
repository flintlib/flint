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

    mask = mpoly_overflow_mask_sp(bits);

    TMP_START;

    pmin = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_set(pmin, poly_exps + N*(len - 1), N);

    if (!mpoly_monomial_is_zero(pmin, N))
    {
        for (i = 0; i < len - 1; i++)
            mpoly_monomial_min(pmin, pmin, poly_exps + N*i, bits, N, mask);
    }

    /* GCC really wants to complain about this one */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    mpoly_unpack_vec_ui(min_fields, pmin, bits, mctx->nfields, 1);
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif

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
    mpoly_monomial_set(pmin, poly_exps + N*(len - 1), N);

    if (!mpoly_monomial_is_zero(pmin, N))
    {
        if (bits <= FLINT_BITS)
        {
            mask = mpoly_overflow_mask_sp(bits);

            for (i = 0; i < len - 1; i++)
                mpoly_monomial_min(pmin, pmin, poly_exps + N*i, bits, N, mask);
        }
        else
        {
            for (i = 0; i < len - 1; i++)
                mpoly_monomial_min_mp(pmin, pmin, poly_exps + N*i, bits, N);
        }
    }

    /* GCC really wants to complain about this one */
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    mpoly_unpack_vec_fmpz(min_fields, pmin, bits, mctx->nfields, 1);
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif

    TMP_END;
}
