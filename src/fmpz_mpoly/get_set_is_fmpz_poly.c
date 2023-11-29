/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mpoly.h"

int fmpz_mpoly_is_fmpz_poly(
    const fmpz_mpoly_t A,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_is_poly(A->exps, A->length, A->bits, var, ctx->minfo);
}

int fmpz_mpoly_get_fmpz_poly(
    fmpz_poly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;

    fmpz_poly_zero(A);

    if (B->length < 1)
        return 1;

    if (Bbits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
        slong off, shift;

        mpoly_gen_offset_shift_sp(&off, &shift, var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = (Bexps[N*i + off] >> shift) & mask;
            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
    else
    {
        slong j, off;
        ulong check, wpf = Bbits/FLINT_BITS;

        off = mpoly_gen_offset_mp(var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = Bexps[N*i + off + 0];
            check = 0;
            for (j = 1; j < wpf; j++)
                check |= Bexps[N*i + off + j];

            if (check != 0 || (slong) k < 0)
                return 0;

            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
}

void _fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz * Bcoeffs,
    slong Blen,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (Bcoeffs[i] != 0);

    fmpz_mpoly_fit_length_reset_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (fmpz_is_zero(Bcoeffs + i))
            continue;

        FLINT_ASSERT(Alen < A->alloc);
        fmpz_set(A->coeffs + Alen, Bcoeffs + i);
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}

void fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    const fmpz_poly_t B,
    slong v,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
    }
    else if (B->length == 1)
    {
        fmpz_mpoly_set_fmpz(A, B->coeffs, ctx);
    }
    else
    {
        bits = mpoly_gen_pow_exp_bits_required(v, B->length - 1, ctx->minfo);
        bits = mpoly_fix_bits(bits, ctx->minfo);
        _fmpz_mpoly_set_fmpz_poly(A, bits, B->coeffs, B->length, v, ctx);
    }
}

void
fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_t res, slong var, const fmpz_poly_t pol, const fmpz_mpoly_ctx_t ctx)
{
    if (ctx->minfo->nvars == 0)
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_set_gen_fmpz_poly: require nvars >= 1");
    }

    fmpz_mpoly_set_fmpz_poly(res, pol, var, ctx);
}

/*
    construct a polynomial with at least Aminbits bits from the coefficients
    Acoeffs[0], ..., Acoeffs[deg]. *** These fmpz's are cleared *** .
*/
void _fmpz_mpoly_set_fmpz_poly_one_var(
    fmpz_mpoly_t A,
    flint_bitcnt_t Aminbits,
    fmpz * Acoeffs,     /* cleared upon return */
    slong Adeg,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, Alen;
    flint_bitcnt_t Abits;

    Abits = mpoly_gen_pow_exp_bits_required(0, Adeg, ctx->minfo);
    Abits = FLINT_MAX(Abits, Aminbits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);
    fmpz_mpoly_fit_length_reset_bits(A, Adeg + 1, Abits, ctx);

    FLINT_ASSERT(Abits <= FLINT_BITS);

    Alen = 0;
    if (ctx->minfo->ord == ORD_LEX)
    {
        FLINT_ASSERT(1 == mpoly_words_per_exp(Abits, ctx->minfo));

        for (i = Adeg; i >= 0; i--)
        {
            if (fmpz_is_zero(Acoeffs + i))
                continue;

            fmpz_swap(A->coeffs + Alen, Acoeffs + i);
            A->exps[Alen] = i;
            Alen++;
            fmpz_clear(Acoeffs + i);
        }
    }
    else if (1 == mpoly_words_per_exp(Abits, ctx->minfo))
    {
        for (i = Adeg; i >= 0; i--)
        {
            if (fmpz_is_zero(Acoeffs + i))
                continue;

            fmpz_swap(A->coeffs + Alen, Acoeffs + i);
            A->exps[Alen] = i + (i << Abits);
            Alen++;
            fmpz_clear(Acoeffs + i);
        }
    }
    else
    {
        FLINT_ASSERT(2 == mpoly_words_per_exp(Abits, ctx->minfo));

        for (i = Adeg; i >= 0; i--)
        {
            if (fmpz_is_zero(Acoeffs + i))
                continue;

            fmpz_swap(A->coeffs + Alen, Acoeffs + i);
            A->exps[2*Alen + 1] = A->exps[2*Alen + 0] = i;
            Alen++;
            fmpz_clear(Acoeffs + i);
        }
    }

    _fmpz_mpoly_set_length(A, Alen, ctx);
}
