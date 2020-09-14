/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx)
{
    slong i, len = b->length;

    if (len <= 2*(ctx->modulus->length - 1))
    {
        fmpz_poly_fit_length(a, len);

        for (i = 0; i < len; i++)
            fmpz_set(a->coeffs + i, b->coeffs + i);

        _fmpz_poly_set_length(a, len);

        fq_reduce(a, ctx);
    }
    else
    {
        fmpz_mod_poly_t q, r;
        fmpz_mod_poly_init(q, ctx->ctxp);
        fmpz_mod_poly_init(r, ctx->ctxp);
        fmpz_mod_poly_divrem(q, r, b, fq_ctx_modulus(ctx), ctx->ctxp);
        fmpz_mod_poly_get_fmpz_poly(a, r, ctx->ctxp);
        fmpz_mod_poly_clear(q, ctx->ctxp);
        fmpz_mod_poly_clear(r, ctx->ctxp);
    }
}
