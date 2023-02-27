/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_make_monic(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    fq_zech_t c;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "Zero polynomial in fq_zech_mpoly_make_monic");
    }

    fq_zech_init(c, ctx->fqctx);
    fq_zech_inv(c, B->coeffs + 0, ctx->fqctx);
    fq_zech_mpoly_scalar_mul_fq_zech(A, B, c, ctx);
    fq_zech_clear(c, ctx->fqctx);
}
