/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

/* exponents of B are not multiprecision */
static void _fq_zech_mpoly_evaluate_one_fq_zech_sp(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_t val,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    slong main_shift, main_off;
    ulong main_exp;
    slong Blen = B->length;
    fq_zech_struct * Bcoeffs = B->coeffs;
    ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong * one;
    fq_zech_struct * Acoeffs;
    ulong * Aexps;
    fq_zech_t pp;
    TMP_INIT;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    fq_zech_init(pp, ctx->fqctx);

    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_offset_shift_sp(one, &main_off, &main_shift,
                                                        var, bits, ctx->minfo);

    fq_zech_mpoly_fit_length_set_bits(A, B->length, bits, ctx);
    Acoeffs = A->coeffs;
    Aexps = A->exps;
    for (i = 0; i < Blen; i++)
    {
        main_exp = (Bexps[N*i + main_off] >> main_shift) & mask;
        fq_zech_pow_ui(pp, val, main_exp, ctx->fqctx);
        fq_zech_mul(Acoeffs + i, Bcoeffs + i, pp, ctx->fqctx);
        mpoly_monomial_msub(Aexps + N*i, Bexps + N*i, main_exp, one, N);
    }

    A->length = Blen;

    fq_zech_mpoly_sort_terms(A, ctx);
    fq_zech_mpoly_combine_like_terms(A, ctx);

    TMP_END;

    return;
}


/* exponents of B are multiprecision */
static void _fq_zech_mpoly_evaluate_one_fq_zech_mp(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_t val,
    const fq_zech_mpoly_ctx_t ctx)
{
    flint_printf("_fq_zech_mpoly_evaluate_one_fq_zech_mp not implemented\n");
    flint_abort();
}

void fq_zech_mpoly_evaluate_one_fq_zech(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_t val,
    const fq_zech_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fq_zech_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        fq_zech_mpoly_t T;
        fq_zech_mpoly_init(T, ctx);
        fq_zech_mpoly_evaluate_one_fq_zech(T, B, var, val, ctx);
        fq_zech_mpoly_swap(A, T, ctx);
        fq_zech_mpoly_clear(T, ctx);
        return;
    }

    if (B->bits <= FLINT_BITS)
    {
        _fq_zech_mpoly_evaluate_one_fq_zech_sp(A, B, var, val, ctx);
    }
    else
    {
        _fq_zech_mpoly_evaluate_one_fq_zech_mp(A, B, var, val, ctx);
    }

    fq_zech_mpoly_assert_canonical(A, ctx);
}

