/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

slong _fmpz_mpoly_derivative(fmpz * coeff1, ulong * exp1,
                       const fmpz * coeff2, const ulong * exp2, slong len2,
     slong var, slong bits, slong N, slong offset, slong shift, ulong * oneexp)
{
    slong i, len1;

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        ulong c = (exp2[N*i + offset] >> shift) & mask;
        if (c != 0)
        {
            mpoly_monomial_sub(exp1 + N*len1, exp2 + N*i, oneexp, N);
            fmpz_mul_ui(coeff1 + len1, coeff2 + i, c);
            len1++;
        }
    }

    return len1;
}

void fmpz_mpoly_derivative(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong bits, N, offset, shift;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;
    bits = poly2->bits;
    fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift,
                                                     var, N, bits, ctx->minfo);

    len1 = _fmpz_mpoly_derivative(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                          var, bits, N, offset, shift, oneexp);

    _fmpz_mpoly_set_length(poly1, len1, ctx);
}
