/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_ctx_init_rand(fq_nmod_mpoly_ctx_t ctx, flint_rand_t state,
                    slong max_nvars, flint_bitcnt_t p_bits_bound, slong deg_bound)
{
    ulong p;
    slong d;
    flint_bitcnt_t p_bits;
    nmod_poly_t poly;

    d = 1 + n_randint(state, deg_bound);

    p_bits = 1 + n_randint(state, p_bits_bound);

    if (p_bits >= FLINT_BITS)
    {
        p = n_randlimb(state);
        p |= UWORD(1) << (FLINT_BITS - 1);
        p &= ~(UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX));
    }
    else
    {
        p = n_randtest_bits(state, p_bits);
    }
    p = n_nextprime(p, 1);

    nmod_poly_init2(poly, p, d + 1);
    nmod_poly_randtest_monic_irreducible(poly, state, d + 1);
    fq_nmod_ctx_init_modulus(ctx->fqctx, poly, "#");
    nmod_poly_clear(poly);

    mpoly_ctx_init_rand(ctx->minfo, state, max_nvars);
}
