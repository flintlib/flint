/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/* exponents of B are not multiprecision */
void _fmpz_mod_mpoly_evaluate_one_fmpz_mod_sp(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_t val,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N, off, shift;
    ulong * cmpmask, * one;
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    fmpz * Acoeffs;
    ulong * Aexps;
    ulong mask, k;
    int need_sort = 0, cmp;
    fmpz_t t;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;

    fmpz_init(t);

    fmpz_mod_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);
    Acoeffs = A->coeffs;
    Aexps = A->exps;

    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_offset_shift_sp(one, &off, &shift, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        k = (Bexps[N*i + off] >> shift) & mask;
        fmpz_mod_pow_ui(t, val, k, ctx->ffinfo);
        fmpz_mod_mul(Acoeffs + Alen, Bcoeffs + i, t, ctx->ffinfo);

        if (fmpz_is_zero(Acoeffs + Alen))
            continue;

        mpoly_monomial_msub(Aexps + N*Alen, Bexps + N*i, k, one, N);
        if (Alen < 1)
        {
            Alen = 1;
            continue;
        }
        cmp = mpoly_monomial_cmp(Aexps + N*(Alen - 1), Aexps + N*Alen, N, cmpmask);
        if (cmp != 0)
        {
            need_sort |= (cmp < 0);
            Alen++;
            continue;
        }
        fmpz_mod_add(Acoeffs + Alen - 1, Acoeffs + Alen - 1, Acoeffs + Alen, ctx->ffinfo);
        Alen -= fmpz_is_zero(Acoeffs + Alen - 1);
    }
    A->length = Alen;

    fmpz_clear(t);

    TMP_END;

    if (need_sort)
    {
        fmpz_mod_mpoly_sort_terms(A, ctx);
        fmpz_mod_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(fmpz_mod_mpoly_is_canonical(A, ctx));
}


/* exponents of B are multiprecision */
static void _fmpz_mod_mpoly_evaluate_one_fmpz_mod_mp(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_t val,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N, off;
    ulong * cmpmask, * one, * tmp;
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    fmpz * Acoeffs;
    ulong * Aexps;
    fmpz_t k, t;
    int need_sort = 0, cmp;
    TMP_INIT;

    FLINT_ASSERT((B->bits % FLINT_BITS) == 0);

    TMP_START;

    fmpz_init(k);
    fmpz_init(t);

    fmpz_mod_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);
    Acoeffs = A->coeffs;
    Aexps = A->exps;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong *) TMP_ALLOC(3*N*sizeof(ulong));
    cmpmask = one + N;
    tmp = cmpmask + N;
    off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(k, Bexps + N*i + off, bits/FLINT_BITS);
        fmpz_mod_pow_fmpz(t, val, k, ctx->ffinfo);
        fmpz_mod_mul(Acoeffs + Alen, Bcoeffs + i, t, ctx->ffinfo);

        if (fmpz_is_zero(Acoeffs + Alen))
            continue;

        mpoly_monomial_mul_fmpz(tmp, one, N, k);
        mpoly_monomial_sub_mp(Aexps + N*Alen, Bexps + N*i, tmp, N);
        if (Alen < 1)
        {
            Alen = 1;
            continue;
        }
        cmp = mpoly_monomial_cmp(Aexps + N*(Alen - 1), Aexps + N*Alen, N, cmpmask);
        if (cmp != 0)
        {
            need_sort |= (cmp < 0);
            Alen++;
            continue;
        }
        fmpz_mod_add(Acoeffs + Alen - 1, Acoeffs + Alen - 1, Acoeffs + Alen, ctx->ffinfo);
        Alen -= fmpz_is_zero(Acoeffs + Alen - 1);
    }
    A->length = Alen;

    fmpz_clear(k);
    fmpz_clear(t);

    TMP_END;

    if (need_sort)
    {
        fmpz_mod_mpoly_sort_terms(A, ctx);
        fmpz_mod_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(fmpz_mod_mpoly_is_canonical(A, ctx));
}


void fmpz_mod_mpoly_evaluate_one_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_t val,
    const fmpz_mod_mpoly_ctx_t ctx)
{
#if FLINT_WANT_ASSERT
    flint_bitcnt_t bits = B->bits;
#endif

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    if (fmpz_mod_is_canonical(val, ctx->ffinfo))
    {
        if (B->bits <= FLINT_BITS)
            _fmpz_mod_mpoly_evaluate_one_fmpz_mod_sp(A, B, var, val, ctx);
        else
            _fmpz_mod_mpoly_evaluate_one_fmpz_mod_mp(A, B, var, val, ctx);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mod_set_fmpz(t, val, ctx->ffinfo);

        if (B->bits <= FLINT_BITS)
            _fmpz_mod_mpoly_evaluate_one_fmpz_mod_sp(A, B, var, t, ctx);
        else
            _fmpz_mod_mpoly_evaluate_one_fmpz_mod_mp(A, B, var, t, ctx);

        fmpz_clear(t);
    }

    /* unwritten rule that output bits = input bits if input != 0 */
    FLINT_ASSERT(A->bits == bits);
}
