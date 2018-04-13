/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_pushback_term_fmpq_fmpz(fmpq_mpoly_t poly,
                const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
{
    mp_bitcnt_t exp_bits;
    fmpz_mpoly_struct * zpoly;
    fmpq_t t;
    slong N;

    fmpq_init(t);

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->zctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->zctx->minfo);
    fmpq_mpoly_fit_bits(poly, exp_bits, ctx);

    zpoly = poly->zpoly;
    N = mpoly_words_per_exp(zpoly->bits, ctx->zctx->minfo);

    fmpz_mpoly_fit_length(zpoly, zpoly->length + 1, ctx->zctx);

    mpoly_set_monomial_pfmpz(zpoly->exps + N*zpoly->length, exp, zpoly->bits,
                                                             ctx->zctx->minfo);

    if (!fmpz_is_one(fmpq_numref(poly->content)))
    {
        _fmpz_vec_scalar_mul_fmpz(zpoly->coeffs, zpoly->coeffs,
                                    zpoly->length, fmpq_numref(poly->content));
        fmpz_one(fmpq_numref(poly->content));
    }

    fmpq_mul_fmpz(t, c, fmpq_denref(poly->content));
    if (!fmpz_is_one(fmpq_denref(t)))
    {
        _fmpz_vec_scalar_mul_fmpz(zpoly->coeffs, zpoly->coeffs,
                                                zpoly->length, fmpq_denref(t));
        fmpz_mul(fmpq_denref(poly->content), fmpq_denref(poly->content),
                                                               fmpq_denref(t));
    }
    fmpz_swap(zpoly->coeffs + zpoly->length, fmpq_numref(t));

    zpoly->length++; /* safe because length is increasing */

    fmpq_clear(t);
}
