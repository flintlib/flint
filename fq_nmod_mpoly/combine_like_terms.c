/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    assuming that the exponents are valid and sorted,
    put the polynomial in canonical form
    i.e.
        2*x^e + 3*x^e -> 5x^e
        2*x^e - 2*x^e -> 0
*/
void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong in, out, N = mpoly_words_per_exp(A->bits, ctx->minfo);

    out = -WORD(1);

    for (in = WORD(0); in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= WORD(0)
            && mpoly_monomial_equal(A->exps + N*out, A->exps + N*in, N))
        {
            fq_nmod_add(A->coeffs + out, A->coeffs + out, A->coeffs + in,
                                                                  ctx-> fqctx);
        }
        else
        {
            if (out < WORD(0) || !fq_nmod_is_zero(A->coeffs + out, ctx->fqctx))
                out++;

            if (out != in)
            {
                mpoly_monomial_set(A->exps + N*out, A->exps + N*in, N);
                fq_nmod_swap(A->coeffs + out, A->coeffs + in, ctx-> fqctx);
            }
        }
    }

    if (out < WORD(0) || !fq_nmod_is_zero(A->coeffs + out, ctx->fqctx))
        out++;

    A->length = out;
}
