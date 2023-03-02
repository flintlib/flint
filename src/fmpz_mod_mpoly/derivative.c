/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


static slong _fmpz_mod_mpoly_derivative(
    fmpz * coeff1, ulong * exp1,
    const fmpz * coeff2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    slong N,
    slong offset,
    slong shift,
    ulong * oneexp,
    const fmpz_mod_ctx_t fctx)
{
    slong i, len1;

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        ulong c = (exp2[N*i + offset] >> shift) & mask;
        if (c == 0)
            continue;
        fmpz_mod_mul_ui(coeff1 + len1, coeff2 + i, c, fctx);
        if (fmpz_is_zero(coeff1 + len1))
            continue;
        mpoly_monomial_sub(exp1 + N*len1, exp2 + N*i, oneexp, N);
        len1++;
    }

    return len1;
}


static slong _fmpz_mod_mpoly_derivative_mp(
    fmpz * coeff1, ulong * exp1,
    const fmpz * coeff2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    slong N,
    slong offset,
    ulong * oneexp,
    const fmpz_mod_ctx_t fctx)
{
    slong i, len1;
    fmpz_t c;
    fmpz_init(c);

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
        if (fmpz_is_zero(c))
            continue;
        fmpz_mod_mul_fmpz(coeff1 + len1, coeff2 + i, c, fctx);
        if (fmpz_is_zero(coeff1 + len1))
            continue;
        mpoly_monomial_sub_mp(exp1 + N*len1, exp2 + N*i, oneexp, N);
        len1++;
    }

    fmpz_clear(c);

    return len1;
}


void fmpz_mod_mpoly_derivative(
    fmpz_mod_mpoly_t poly1,
    const fmpz_mod_mpoly_t poly2,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = poly2->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong offset, shift;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;

    fmpz_mod_mpoly_fit_length_reset_bits(poly1, poly2->length, bits, ctx);

    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift,
                                                        var, bits, ctx->minfo);

        len1 = _fmpz_mod_mpoly_derivative(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                  bits, N, offset, shift, oneexp, ctx->ffinfo);
    }
    else
    {
        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, ctx->minfo);

        len1 = _fmpz_mod_mpoly_derivative_mp(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                  bits, N, offset,        oneexp, ctx->ffinfo);
    }

    _fmpz_mod_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
