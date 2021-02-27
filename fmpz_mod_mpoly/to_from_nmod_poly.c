/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/*
    set A(x_var^Bstride[var]) to B/xbar^Bshifts
    it is asserted that the conversion is correct
*/
void _fmpz_mod_mpoly_to_fmpz_mod_poly_deflate(
    fmpz_mod_poly_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const ulong * Bshift,
    const ulong * Bstride,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    ulong mask;
    slong i, shift, off, N;
    slong len = B->length;
    const fmpz * coeff = B->coeffs;
    ulong * exp = B->exps;
    ulong var_shift, var_stride;
    flint_bitcnt_t bits = B->bits;

    FLINT_ASSERT(len > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    fmpz_mod_poly_zero(A, ctx->ffinfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    var_shift = Bshift[var];
    var_stride = Bstride[var];
    for (i = 0; i < len; i++)
    {
        ulong k = (exp[N*i + off] >> shift) & mask;
        FLINT_ASSERT(k >= var_shift);
        k -= var_shift;
        if (k != 0)
            k /= var_stride;
        fmpz_mod_poly_set_coeff_fmpz(A, k, coeff + i, ctx->ffinfo);
    }

#if FLINT_WANT_ASSERT
    for (i = 0; i < len; i++)
    {
        slong v;
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            ulong k;
            mpoly_gen_offset_shift_sp(&off, &shift, v, bits, ctx->minfo);
            k = (exp[N*i + off] >> shift) & mask;
            FLINT_ASSERT(   (v == var && k >= Bshift[v])
                         || (v != var && k == Bshift[v]));
        }
    }
#endif
}

/*
    set A to B(x_var^Astride[var])*xbar^Ashift
    A must be packed into bits = Abits
*/
void _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_poly_t B,
    slong var,
    const ulong * Ashift,
    const ulong * Astride,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    ulong * shiftexp;
    ulong * strideexp;
    slong Bdeg = fmpz_mod_poly_degree(B, ctx->ffinfo);
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(!fmpz_mod_poly_is_zero(B, ctx->ffinfo));

    /* must have at least space for the highest exponent of var */
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Ashift[var] + Bdeg*Astride[var]) <= Abits);

    N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    shiftexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ui(shiftexp, Ashift, Abits, ctx->minfo);
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);
    mpoly_monomial_mul_ui(strideexp, strideexp, N, Astride[var]);

    fmpz_mod_mpoly_fit_length_reset_bits(A, 0, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (k = Bdeg; k >= 0; k--)
    {
        _fmpz_mod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, N, Alen + 1);
        fmpz_mod_poly_get_coeff_fmpz(Acoeff + Alen, B, k, ctx->ffinfo);
        if (fmpz_is_zero(Acoeff + Alen))
            continue;
        mpoly_monomial_madd(Aexp + N*Alen, shiftexp, k, strideexp, N);
        Alen++;
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    _fmpz_mod_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}
