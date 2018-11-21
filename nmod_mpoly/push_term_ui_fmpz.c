/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    emplaceterm assumes that c is valid modulo ctx->ffinfo->mod.n
*/

void _nmod_mpoly_emplacebackterm_ui_ffmpz(nmod_mpoly_t A,
                   mp_limb_t c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    mp_bitcnt_t exp_bits;
    FLINT_ASSERT(c < ctx->ffinfo->mod.n);

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, old_length + 1, ctx);
    A->length = old_length + 1;
    A->coeffs[old_length] = c;
    mpoly_set_monomial_ffmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}

void _nmod_mpoly_emplacebackterm_ui_pfmpz(nmod_mpoly_t A,
                   mp_limb_t c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    mp_bitcnt_t exp_bits;
    FLINT_ASSERT(c < ctx->ffinfo->mod.n);

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, old_length + 1, ctx);
    A->length = old_length + 1;
    A->coeffs[old_length] = c;
    mpoly_set_monomial_pfmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}


void nmod_mpoly_push_term_ui_fmpz(nmod_mpoly_t A, ulong c,
                                fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
{
    if (c >= ctx->ffinfo->mod.n)
    {
        NMOD_RED(c, c, ctx->ffinfo->mod);
    }
    _nmod_mpoly_emplacebackterm_ui_pfmpz(A, c, exp, ctx);    
}
