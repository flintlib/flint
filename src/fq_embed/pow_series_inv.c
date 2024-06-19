/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_embed.h"

void
fq_modulus_pow_series_inv(fmpz_mod_poly_t res, const fq_ctx_t ctx, slong trunc)
{
    fmpz_mod_poly_reverse(res, fq_ctx_modulus(ctx), fq_ctx_degree(ctx) + 1, ctx->ctxp);
    fmpz_mod_poly_inv_series(res, res, trunc, ctx->ctxp);
}
