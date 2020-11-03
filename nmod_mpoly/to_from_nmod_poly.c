/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    set A(x_var^Bstride[var]) to B/xbar^Bshifts
    it is asserted that the conversion is correct
*/
void _nmod_mpoly_to_nmod_poly_deflate(
    nmod_poly_t A,
    const nmod_mpoly_t B,
    slong var,
    const ulong * Bshift,
    const ulong * Bstride,
    const nmod_mpoly_ctx_t ctx)
{
    ulong mask;
    slong i, shift, off, N;
    slong len = B->length;
    mp_limb_t * coeff = B->coeffs;
    ulong * exp = B->exps;
    ulong var_shift, var_stride;
    flint_bitcnt_t bits = B->bits;

    FLINT_ASSERT(len > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    nmod_poly_zero(A);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    var_shift = Bshift[var];
    var_stride = Bstride[var];
    for (i = 0; i < len; i++)
    {
        ulong k = (exp[N*i + off] >> shift) & mask;
        FLINT_ASSERT(k >= var_shift);
        k -= var_shift;
        if (k != 0)
        {
            k /= var_stride;
        }
        nmod_poly_set_coeff_ui(A, k, coeff[i]);
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
void _nmod_mpoly_from_nmod_poly_inflate(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_poly_t B,
    slong var,
    const ulong * Ashift,
    const ulong * Astride,
    const nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong k;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    ulong * shiftexp;
    ulong * strideexp;
    slong Bdeg = nmod_poly_degree(B);
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(!nmod_poly_is_zero(B));

    /* must have at least space for the highest exponent of var */
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Ashift[var] + Bdeg*Astride[var]) <= Abits);

    N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    shiftexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ui(shiftexp, Ashift, Abits, ctx->minfo);
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);
    mpoly_monomial_mul_ui(strideexp, strideexp, N, Astride[var]);

    nmod_mpoly_fit_length_reset_bits(A, 0, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (k = Bdeg; k >= 0; k--)
    {
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, N, Alen + 1);
        Acoeff[Alen] = nmod_poly_get_coeff_ui(B, k);
        if (Acoeff[Alen] == 0)
            continue;
        mpoly_monomial_madd(Aexp + N*Alen, shiftexp, k, strideexp, N);
        Alen++;
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    _nmod_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}
