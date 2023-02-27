/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_randtest_bits(
    fq_nmod_mpoly_t A,
    flint_rand_t state,
    slong length,
    flint_bitcnt_t exp_bits,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j, nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = mpoly_fix_bits(FLINT_MAX(exp_bits, 1), ctx->minfo);
    fmpz * exp;
    TMP_INIT;

    TMP_START;

    exp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exp + j);

    fq_nmod_mpoly_fit_length_reset_bits(A, length, bits, ctx);
    A->length = 0;
    for (i = 0; i < length; i++)
    {
        mpoly_monomial_randbits_fmpz(exp, state, exp_bits, ctx->minfo);
        _fq_nmod_mpoly_push_exp_ffmpz(A, exp, ctx);
        n_fq_randtest_not_zero(A->coeffs + d*(A->length - 1), state, ctx->fqctx);
    }
    fq_nmod_mpoly_sort_terms(A, ctx);
    fq_nmod_mpoly_combine_like_terms(A, ctx);

    for (j = 0; j < nvars; j++)
        fmpz_clear(exp + j);

    TMP_END;
}
