/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

slong _fmpz_mpoly_integral(fmpz_t s, fmpz * coeff1, ulong * exp1,
                         const fmpz * coeff2, const ulong * exp2, slong len2,
                                 slong var, slong bits, const mpoly_ctx_t mctx)
{
    slong i, N;
    fmpz_t d, g;
    slong offset;
    slong shift;
    ulong * oneexp;
    TMP_INIT;

    TMP_START;
    fmpz_init(d);
    fmpz_init(g);
    fmpz_set_ui(s, WORD(1));

    N = mpoly_words_per_exp(bits, mctx);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (bits <= FLINT_BITS)
    {
        ulong c, mask;

        mask = (-UWORD(1)) >> (FLINT_BITS - bits);        
        mpoly_gen_monomial_offset_shift_sp(oneexp, &offset, &shift,
                                                              var, bits, mctx);

        /* scan once to find required denominator */
        for (i = 0; i < len2; i++)
        {
            c = (exp2[N*i + offset] >> shift) & mask;
            fmpz_set_ui(d, c + 1);
            fmpz_gcd(g, d, coeff2 + i);
            fmpz_divexact(d, d, g);
            fmpz_lcm(s, s, d);
        }

        /* then scan again to compute the terms */
        /* x^c -> x^(c+1)/(c+1) */
        for (i = 0; i < len2; i++)
        {
            c = (exp2[N*i + offset] >> shift) & mask;
            mpoly_monomial_add(exp1 + N*i, exp2 + N*i, oneexp, N);
            fmpz_set_ui(d, c + 1);
            fmpz_mul(g, coeff2 + i, s);
            fmpz_mul(coeff1 + i, coeff2 + i, g);
            fmpz_divexact(coeff1 + i, g, d);
        }

    } else
    {
        fmpz_t c;
        fmpz_init(c);

        offset = mpoly_gen_monomial_offset_mp(oneexp, var, bits, mctx);

        /* scan once to find required denominator */
        for (i = 0; i < len2; i++)
        {
            fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
            fmpz_add_ui(d, c, 1);
            fmpz_gcd(g, d, coeff2 + i);
            fmpz_divexact(d, d, g);
            fmpz_lcm(s, s, d);
        }

        /* then scan again to compute the terms */
        /* x^c -> x^(c+1)/(c+1) */
        for (i = 0; i < len2; i++)
        {
            fmpz_set_ui_array(c, exp2 + N*i + offset, bits/FLINT_BITS);
            mpoly_monomial_add_mp(exp1 + N*i, exp2 + N*i, oneexp, N);
            fmpz_add_ui(d, c, 1);
            fmpz_mul(g, coeff2 + i, s);
            fmpz_mul(coeff1 + i, coeff2 + i, g);
            fmpz_divexact(coeff1 + i, g, d);
        }

        fmpz_clear(c);
    }

    fmpz_clear(g);
    fmpz_clear(d);
    TMP_END;
    return len2;
}

void fmpz_mpoly_integral(fmpz_mpoly_t poly1, fmpz_t scale,
               const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, len1, exp_bits;
    ulong * exp2 = poly2->exps;
    fmpz * gen_fields, * max_fields;
    int free2 = 0;
    TMP_INIT;

    TMP_START;

    /* compute bits required to represent result */
    gen_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(gen_fields + i);
        fmpz_init(max_fields + i);
    }
    mpoly_gen_fields_fmpz(gen_fields, var, ctx->minfo);
    mpoly_max_fields_fmpz(max_fields, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    _fmpz_vec_add(max_fields, max_fields, gen_fields, ctx->minfo->nfields);

    exp_bits = _fmpz_vec_max_bits(max_fields, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(gen_fields + i);
        fmpz_clear(max_fields + i);
    }


    /* ensure input exponents are packed into same sized fields as output */
    if (exp_bits > poly2->bits)
    {
        slong N = mpoly_words_per_exp(exp_bits, ctx->minfo);
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
    }

    /* deal with aliasing and do integration */
    if (poly1 == poly2)
    {
        fmpz_mpoly_t temp;

        fmpz_mpoly_init2(temp, poly2->length, ctx);
        fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;

        len1 = _fmpz_mpoly_integral(scale, temp->coeffs, temp->exps,
                                          poly2->coeffs, exp2, poly2->length,
                                                    var, exp_bits, ctx->minfo);
        _fmpz_mpoly_set_length(temp, len1, ctx);

        fmpz_mpoly_swap(temp, poly1, ctx);
        fmpz_mpoly_clear(temp, ctx);
    } else
    {
        fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len1 = _fmpz_mpoly_integral(scale, poly1->coeffs, poly1->exps,
                                           poly2->coeffs, exp2, poly2->length,
                                                    var, exp_bits, ctx->minfo);
        _fmpz_mpoly_set_length(poly1, len1, ctx);
    }

    if (free2)
    {
        flint_free(exp2);
    }

    TMP_END;
}
