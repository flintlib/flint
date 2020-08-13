/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

slong _fq_zech_mpoly_derivative(fq_zech_struct * Acoeff,       ulong * Aexp,
                 const fq_zech_struct * Bcoeff, const ulong * Bexp, slong Blen,
          flint_bitcnt_t bits, slong N, slong offset, slong shift, ulong * oneexp,
                                                     const fq_zech_ctx_t fqctx)
{
    slong i, Alen;

    /* x^c -> c*x^(c-1) */
    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        mp_limb_t cr;
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        ulong c = (Bexp[N*i + offset] >> shift) & mask;
        if (c == 0)
            continue;
        NMOD_RED(cr, c, fq_zech_ctx_mod(fqctx));
        if (cr == 0)
            continue;
        fq_zech_mul_ui(Acoeff + Alen, Bcoeff + i, cr, fqctx);
        mpoly_monomial_sub(Aexp + N*Alen, Bexp + N*i, oneexp, N);
        Alen++;
    }

    return Alen;
}


slong _fq_zech_mpoly_derivative_mp(fq_zech_struct * Acoeff,       ulong * Aexp,
                 const fq_zech_struct * Bcoeff, const ulong * Bexp, slong Blen,
          flint_bitcnt_t bits, slong N, slong offset,              ulong * oneexp,
                                                     const fq_zech_ctx_t fqctx)
{
    slong i, Alen;
    fmpz_t c;
    fmpz_init(c);

    /* x^c -> c*x^(c-1) */
    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        mp_limb_t cr;
        fmpz_set_ui_array(c, Bexp + N*i + offset, bits/FLINT_BITS);
        if (fmpz_is_zero(c))
            continue;
        cr = fmpz_fdiv_ui(c, fq_zech_ctx_mod(fqctx).n);
        if (cr == 0)
            continue;
        fq_zech_mul_ui(Acoeff + Alen, Bcoeff + i, cr, fqctx);
        mpoly_monomial_sub_mp(Aexp + N*Alen, Bexp + N*i, oneexp, N);
        Alen++;
    }

    fmpz_clear(c);

    return Alen;
}


void fq_zech_mpoly_derivative(fq_zech_mpoly_t poly1, const fq_zech_mpoly_t poly2,
                                      slong var, const fq_zech_mpoly_ctx_t ctx)
{
    slong bits, N, offset, shift;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;
    bits = poly2->bits;

    fq_zech_mpoly_fit_length(poly1, poly2->length, ctx);
    fq_zech_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift,
                                                        var, bits, ctx->minfo);

        len1 = _fq_zech_mpoly_derivative(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                   bits, N, offset, shift, oneexp, ctx->fqctx);
    }
    else
    {
        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, ctx->minfo);

        len1 = _fq_zech_mpoly_derivative_mp(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                                   bits, N, offset,        oneexp, ctx->fqctx);        
    }

    _fq_zech_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
