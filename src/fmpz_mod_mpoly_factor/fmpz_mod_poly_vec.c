/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


slong _fmpz_mod_poly_vec_max_degree(const fmpz_mod_poly_struct * A,
                                         slong Alen, const fmpz_mod_ctx_t ctx)
{
    slong i, len = 0;
    for (i = 0; i < Alen; i++)
        len = FLINT_MAX(len, fmpz_mod_poly_length(A + i, ctx));
    return len - 1;
}

void _fmpz_mod_poly_vec_content(
    fmpz_mod_poly_t g,
    const fmpz_mod_poly_struct * A,
    slong Alen,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_poly_zero(g, ctx);

    for (i = 0; i < Alen; i++)
    {
        fmpz_mod_poly_gcd(g, g, A + i, ctx);

        if (fmpz_mod_poly_is_one(g, ctx))
            break;
    }
}

void _fmpz_mod_poly_vec_remove_content(
    fmpz_mod_poly_t g,
    fmpz_mod_poly_struct * A,
    slong Alen,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t r;

    _fmpz_mod_poly_vec_content(g, A, Alen, ctx);

    if (fmpz_mod_poly_is_one(g, ctx))
        return;

    fmpz_mod_poly_init(r, ctx);

    for (i = 0; i < Alen; i++)
        fmpz_mod_poly_divrem(A + i, r, A + i, g, ctx);

    fmpz_mod_poly_clear(r, ctx);
}

void _fmpz_mod_poly_vec_mul_poly(
    fmpz_mod_poly_struct * A,
    slong Alen,
    const fmpz_mod_poly_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (fmpz_mod_poly_is_one(g, ctx))
        return;

    for (i = 0; i < Alen; i++)
        fmpz_mod_poly_mul(A + i, A + i, g, ctx);
}

void _fmpz_mod_poly_vec_divexact_poly(
    fmpz_mod_poly_struct * A,
    slong Alen,
    const fmpz_mod_poly_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t t;

    if (fmpz_mod_poly_is_one(g, ctx))
        return;

    fmpz_mod_poly_init(t, ctx);

    for (i = 0; i < Alen; i++)
    {
        fmpz_mod_poly_divrem(A + i, t, A + i, g, ctx);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(t, ctx));
    }

    fmpz_mod_poly_clear(t, ctx);
}

void _fmpz_mod_poly_vec_mul_fmpz_mod(
    fmpz_mod_poly_struct * A,
    slong Alen,
    const fmpz_t g,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (fmpz_is_one(g))
        return;

    for (i = 0; i < Alen; i++)
        fmpz_mod_poly_scalar_mul_fmpz(A + i, A + i, g, ctx);
}

