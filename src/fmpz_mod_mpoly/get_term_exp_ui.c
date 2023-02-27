/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mod_mpoly_t A, 
                                       slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_get_term_exp_ui: index out of range");

    mpoly_get_monomial_ui(exp, A->exps + N*i, A->bits, ctx->minfo);
}
