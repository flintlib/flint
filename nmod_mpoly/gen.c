/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_gen(nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
{
    slong j;
    ulong * mon;
    mp_bitcnt_t bits;
    TMP_INIT;

    TMP_START;

    mon = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    for (j = 0; j < ctx->minfo->nvars; j++)
       mon[j] = (j == var);

    nmod_mpoly_fit_length(A, WORD(1), ctx);
    A->coeffs[0] = UWORD(1);
    _nmod_mpoly_set_length(A, WORD(1), ctx);

    bits = mpoly_exp_bits_required_ui(mon, ctx->minfo);
    nmod_mpoly_fit_bits(A, bits, ctx);
    mpoly_set_monomial_ui(A->exps, mon, A->bits, ctx->minfo);

    TMP_END;
}
