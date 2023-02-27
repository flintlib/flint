/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_set_term_exp_ui(fmpz_mod_mpoly_t A, 
                   slong i, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t exp_bits;

    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_set_term_exp_ui: index is out of range");

    exp_bits = mpoly_exp_bits_required_ui(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mod_mpoly_fit_length_fit_bits(A, A->length, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_ui(A->exps + N*i, exp, A->bits, ctx->minfo);
}
