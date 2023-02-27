/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"

void fmpz_mod_pow_cache_start(
    const fmpz_t b,
    fmpz_mod_poly_t c,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(c, 2, ctx);
    c->length = 2;
    fmpz_one(c->coeffs + 0);
    fmpz_set(c->coeffs + 1, b);
}

/* a = b*c^e */
void fmpz_mod_pow_cache_mulpow_ui(
    fmpz_t a,
    const fmpz_t b,
    ulong e,
    fmpz_mod_poly_t c,
    const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(c->length > 1);

    if (e < c->length)
    {
        fmpz_mod_mul(a, b, c->coeffs + e, ctx);
        return;
    }

    if (e > 50)
    {
        fmpz_mod_poly_fit_length(c, c->length + 1, ctx);
        fmpz_mod_pow_ui(c->coeffs + c->length, c->coeffs + 1, e, ctx);
        fmpz_mod_mul(a, b, c->coeffs + c->length, ctx);
        return;
    }

    fmpz_mod_poly_fit_length(c, e + 1, ctx);

    while (e >= c->length)
    {
        fmpz_mod_mul(c->coeffs + c->length,
                     c->coeffs + c->length - 1, c->coeffs + 1, ctx);
        c->length++;
    }

    fmpz_mod_mul(a, b, c->coeffs + e, ctx);
    return;
}

