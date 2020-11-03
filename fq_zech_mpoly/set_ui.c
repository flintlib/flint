/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_set_ui(fq_zech_mpoly_t A, ulong c,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);

    fq_zech_mpoly_fit_length_reset_bits(A, 1, bits, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    fq_zech_set_ui(A->coeffs + 0, c, ctx->fqctx);
    mpoly_monomial_zero(A->exps, N);
    A->length = !fq_zech_is_zero(A->coeffs + 0, ctx->fqctx);
}
