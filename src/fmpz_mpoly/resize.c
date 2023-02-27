/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_resize(fmpz_mpoly_t A, slong new_length,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong old_length = A->length;
    slong N;

    new_length = FLINT_MAX(WORD(0), new_length);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (new_length < old_length)
    {
        /* just zero out the extra unused fmpz's past the new end */
        _fmpz_vec_zero(A->coeffs + new_length, old_length - new_length);
    }
    else if (new_length > old_length)
    {
        if (new_length > A->alloc)
            fmpz_mpoly_realloc(A, FLINT_MAX(new_length, 2*A->alloc), ctx);

        /* must zero out the new coeffs/exps past the old end */
        flint_mpn_zero(A->exps + N*old_length, N*(new_length - old_length));
        _fmpz_vec_zero(A->coeffs + old_length, new_length - old_length);
    }

    A->length = new_length;
}
