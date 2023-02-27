/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

slong _fmpz_mpoly_derivative(fmpz * coeff1, ulong * exp1,
                       const fmpz * coeff2, const ulong * exp2, slong len2,
          flint_bitcnt_t bits, slong N, slong offset, slong shift, ulong * oneexp)
{
    slong i, len1;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
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


slong _fmpz_mpoly_derivative_mp(fmpz * coeff1, ulong * exp1,
                       const fmpz * coeff2, const ulong * exp2, slong len2,
          flint_bitcnt_t bits, slong N, slong offset,              ulong * oneexp)
{
    slong i, len1;
    fmpz_t c;
    fmpz_init(c);

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
        if (!fmpz_is_zero(c))
        {
            mpoly_monomial_sub_mp(exp1 + N*len1, exp2 + N*i, oneexp, N);
            fmpz_mul(coeff1 + len1, coeff2 + i, c);
            len1++;
        }
    }

    fmpz_clear(c);

    return len1;
}


void fmpz_mpoly_derivative(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong N, offset, shift;
    flint_bitcnt_t bits;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;
    bits = poly2->bits;

    if (poly1 != poly2)
        fmpz_mpoly_fit_length_reset_bits(poly1, poly2->length, bits, ctx);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS) {
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift,
                                                        var, bits, ctx->minfo);

        len1 = _fmpz_mpoly_derivative(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                               bits, N, offset, shift, oneexp);
    } else
    {
        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, ctx->minfo);

        len1 = _fmpz_mpoly_derivative_mp(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                               bits, N, offset,        oneexp);        
    }

    _fmpz_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
