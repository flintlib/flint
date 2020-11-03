/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    set A(x_var^Bstride[var]) to B/xbar^Bshifts
    it is asserted that the conversion is correct
*/
void _fq_nmod_mpoly_to_fq_nmod_poly_deflate(
    fq_nmod_poly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const ulong * Bshift,
    const ulong * Bstride,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    ulong mask;
    slong i, shift, off, N;
    slong len = B->length;
    mp_limb_t * coeff = B->coeffs;
    ulong * exp = B->exps;
    ulong var_shift, var_stride;
    flint_bitcnt_t bits = B->bits;
    fq_nmod_t cc;

    FLINT_ASSERT(len > 0);
    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    fq_nmod_init(cc, ctx->fqctx);

    fq_nmod_poly_zero(A, ctx->fqctx);
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
        n_fq_get_fq_nmod(cc, coeff + d*i, ctx->fqctx);
        fq_nmod_poly_set_coeff(A, k, cc, ctx->fqctx);
    }

    fq_nmod_clear(cc, ctx->fqctx);

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
void _fq_nmod_mpoly_from_fq_nmod_poly_inflate(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_poly_t B,
    slong var,
    const ulong * Ashift,
    const ulong * Astride,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong k, Alen;
    ulong * shiftexp, * strideexp;
    slong Bdeg = fq_nmod_poly_degree(B, ctx->fqctx);
    TMP_INIT;

    TMP_START;

    /* must have at least space for the highest exponent of var */
    FLINT_ASSERT(Bdeg >= 0);
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Ashift[var] + Bdeg*Astride[var]) <= Abits);
    
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    shiftexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ui(shiftexp, Ashift, Abits, ctx->minfo);
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);
    mpoly_monomial_mul_ui(strideexp, strideexp, N, Astride[var]);

    fq_nmod_mpoly_fit_length_reset_bits(A, Bdeg + 1, Abits, ctx);
    Alen = 0;
    for (k = Bdeg; k >= 0; k--)
    {
        n_fq_set_fq_nmod(A->coeffs + d*Alen, B->coeffs + k, ctx->fqctx);
        if (!_n_fq_is_zero(A->coeffs + d*Alen, d))
        {
            mpoly_monomial_madd(A->exps + N*Alen, shiftexp, k, strideexp, N);
            Alen++;
        }
    }
    _fq_nmod_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}
