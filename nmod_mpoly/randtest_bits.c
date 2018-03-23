/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_randtest_bits(nmod_mpoly_t poly, flint_rand_t state,
                slong length, mp_bitcnt_t exp_bits, const nmod_mpoly_ctx_t ctx)
{
    slong i, j, nvars = ctx->minfo->nvars;
    mp_limb_t coeff;
    fmpz * exp;
    TMP_INIT;

    TMP_START;

    exp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exp + j);

    nmod_mpoly_zero(poly, ctx);
    for (i = 0; i < length; i++)
    {
        mpoly_monomial_randbits_fmpz(exp, state, exp_bits, ctx->minfo);
        coeff = n_randint(state, ctx->ffinfo->mod.n);
        _nmod_mpoly_set_term_ui_fmpz(poly, coeff, exp, ctx);
    }

    for (j = 0; j < nvars; j++)
        fmpz_clear(exp + j);
}
