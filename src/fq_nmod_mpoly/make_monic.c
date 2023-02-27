/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_make_monic(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    mp_limb_t * c;
    TMP_INIT;

    if (B->length < 1)
    {
        flint_throw(FLINT_ERROR, "fq_nmod_mpoly_make_monic: input is zero");
    }

    TMP_START;
    c = (mp_limb_t *) TMP_ALLOC((1 + N_FQ_INV_ITCH)*d*sizeof(mp_limb_t));

    _n_fq_inv(c, B->coeffs + d*0, ctx->fqctx, c + d);
    fq_nmod_mpoly_scalar_mul_n_fq(A, B, c, ctx);

    TMP_END;
}
