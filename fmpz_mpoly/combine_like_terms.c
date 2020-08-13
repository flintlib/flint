/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    assuming that the exponents are valid and sorted,
    put the polynomial in canonical form
    i.e.
        2*x^e + 3*x^e -> 5x^e
        2*x^e - 2*x^e -> 0
*/
void fmpz_mpoly_combine_like_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong in, out, N = mpoly_words_per_exp(A->bits, ctx->minfo);

    out = -WORD(1);

    for (in = WORD(0); in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= WORD(0) &&
                     mpoly_monomial_equal(A->exps + N*out, A->exps + N*in, N))
        {
            fmpz_add(A->coeffs + out, A->coeffs + out, A->coeffs + in);
        } else
        {
            if (out < WORD(0) || !fmpz_is_zero(A->coeffs + out))
                out++;

            if (out != in)
            {
                mpoly_monomial_set(A->exps + N*out, A->exps + N*in, N);
                fmpz_swap(A->coeffs + out, A->coeffs + in);
            }
        }
    }

    if (out < WORD(0) || !fmpz_is_zero(A->coeffs + out))
        out++;

    _fmpz_mpoly_set_length(A, out, ctx);
}
