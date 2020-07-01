/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq.h"

void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx)
{
    slong i, len = b->length;

    FLINT_ASSERT(fmpz_equal(&b->p, &ctx->p));

    fmpz_poly_fit_length(a, len);

    for (i = 0; i < len; i++)
        fmpz_set(a->coeffs + i, b->coeffs + i);

    _fmpz_poly_set_length(a, len);

    fq_reduce(a, ctx);
}
