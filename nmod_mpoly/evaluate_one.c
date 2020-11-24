/*
    Copyright (C) 2018-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "n_poly.h"

/* exponents of B are not multiprecision */
void _nmod_mpoly_evaluate_one_ui_sp(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    slong var,
    ulong val,
    const nmod_mpoly_ctx_t ctx,
    n_poly_stack_t St)
{
    slong i, N, off, shift;
    ulong * cmpmask, * one;
    slong Blen = B->length;
    const mp_limb_t * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    mp_limb_t * Acoeffs;
    ulong * Aexps;
    ulong mask, k;
    int need_sort = 0, cmp;
    n_poly_struct * cache[3];
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;

    n_poly_stack_fit_request(St, 3);
    cache[0] = n_poly_stack_take_top(St);
    cache[1] = n_poly_stack_take_top(St);
    cache[2] = n_poly_stack_take_top(St);
    nmod_pow_cache_start(val, cache[0], cache[1], cache[2]);

    nmod_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);
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
        Acoeffs[Alen] = nmod_pow_cache_mulpow_ui(Bcoeffs[i], k, cache[0],
                                                 cache[1], cache[2], ctx->mod);
        if (Acoeffs[Alen] == 0)
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
        Acoeffs[Alen - 1] = nmod_add(Acoeffs[Alen - 1], Acoeffs[Alen], ctx->mod);
        Alen -= (Acoeffs[Alen - 1] == 0);
    }
    A->length = Alen;

    n_poly_stack_give_back(St, 3);

    TMP_END;

    if (need_sort)
    {
        nmod_mpoly_sort_terms(A, ctx);
        nmod_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(nmod_mpoly_is_canonical(A, ctx));
}


/* exponents of B are multiprecision */
static void _nmod_mpoly_evaluate_one_ui_mp(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    slong var,
    ulong val,
    const nmod_mpoly_ctx_t ctx,
    n_poly_stack_t St)
{
    slong i, N, off;
    ulong * cmpmask, * one, * tmp;
    slong Blen = B->length;
    const mp_limb_t * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    mp_limb_t * Acoeffs;
    ulong * Aexps;
    fmpz_t k;
    int need_sort = 0, cmp;
    n_poly_struct * cache[3];
    TMP_INIT;

    FLINT_ASSERT((B->bits % FLINT_BITS) == 0);

    TMP_START;

    fmpz_init(k);

    n_poly_stack_fit_request(St, 3);
    cache[0] = n_poly_stack_take_top(St);
    cache[1] = n_poly_stack_take_top(St);
    cache[2] = n_poly_stack_take_top(St);
    nmod_pow_cache_start(val, cache[0], cache[1], cache[2]);

    nmod_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);
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
        Acoeffs[Alen] = nmod_pow_cache_mulpow_fmpz(Bcoeffs[i], k, cache[0],
                                                 cache[1], cache[2], ctx->mod);
        if (Acoeffs[Alen] == 0)
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
        Acoeffs[Alen - 1] = nmod_add(Acoeffs[Alen - 1], Acoeffs[Alen], ctx->mod);
        Alen -= (Acoeffs[Alen - 1] == 0);
    }
    A->length = Alen;

    n_poly_stack_give_back(St, 3);
    fmpz_clear(k);

    TMP_END;

    if (need_sort)
    {
        nmod_mpoly_sort_terms(A, ctx);
        nmod_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(nmod_mpoly_is_canonical(A, ctx));
}


void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                              slong var, ulong val, const nmod_mpoly_ctx_t ctx)
{
    n_poly_stack_t St;

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    if (val >= ctx->mod.n)
        NMOD_RED(val, val, ctx->mod);

    n_poly_stack_init(St);

    if (B->bits <= FLINT_BITS)
        _nmod_mpoly_evaluate_one_ui_sp(A, B, var, val, ctx, St);
    else
        _nmod_mpoly_evaluate_one_ui_mp(A, B, var, val, ctx, St);

    n_poly_stack_clear(St);
}
