/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

void fq_get_fmpz_mod_poly(fmpz_mod_poly_t a, const fq_t b, const fq_ctx_t ctx)
{
    slong i, len = b->length;

    fmpz_mod_poly_fit_length(a, len, ctx->ctxp);

    for (i = 0; i < len; i++)
        fmpz_set(a->coeffs + i, b->coeffs + i);

    _fmpz_mod_poly_set_length(a, len);
    _fmpz_mod_poly_normalise(a);
}
