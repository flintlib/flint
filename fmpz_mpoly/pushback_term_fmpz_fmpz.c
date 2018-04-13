/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


void fmpz_mpoly_pushback_term_fmpz_fmpz(fmpz_mpoly_t poly,
                 const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    mp_bitcnt_t exp_bits;
    slong N;

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    fmpz_mpoly_fit_length(poly, poly->length + 1, ctx);
    fmpz_set(poly->coeffs + poly->length, c);
    mpoly_set_monomial_pfmpz(poly->exps + N*poly->length, exp, poly->bits, ctx->minfo);
    poly->length++; /* safe because length is increasing */
}
