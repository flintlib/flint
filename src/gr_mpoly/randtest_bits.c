/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_mpoly.h"

int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state,
                 slong length, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong i, j, nvars = mctx->nvars;
    fmpz * exp;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;

    exp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exp + j);

    status = gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length_reset_bits(A, 0, MPOLY_MIN_BITS, mctx, cctx);

    for (i = 0; i < length; i++)
    {
        mpoly_monomial_randbits_fmpz(exp, state, exp_bits, mctx);
        _gr_mpoly_push_exp_fmpz(A, exp, mctx, cctx);
        status |= gr_randtest(GR_ENTRY(A->coeffs, A->length - 1, cctx->sizeof_elem), state, cctx);
    }

    gr_mpoly_sort_terms(A, mctx, cctx);
    status |= gr_mpoly_combine_like_terms(A, mctx, cctx);

    for (j = 0; j < nvars; j++)
        fmpz_clear(exp + j);

    TMP_END;

    return status;
}
