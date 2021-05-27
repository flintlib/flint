/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

void _n_poly_vec_mod_mul_poly(
    n_poly_struct * A,
    slong Alen,
    const n_poly_t g,
    const nmod_t ctx)
{
    slong i;

    if (n_poly_is_one(g))
        return;

    for (i = 0; i < Alen; i++)
        n_poly_mod_mul(A + i, A + i, g, ctx);
}

void _n_poly_vec_mod_content(
    n_poly_t g,
    const n_poly_struct * A,
    slong Alen,
    nmod_t ctx)
{
    slong i;

    n_poly_zero(g);

    for (i = 0; i < Alen; i++)
    {
        n_poly_mod_gcd(g, g, A + i, ctx);
        if (n_poly_is_one(g))
            break;
    }
}

void _n_poly_vec_mod_remove_content(
    n_poly_t g,
    n_poly_struct * A,
    slong Alen,
    nmod_t ctx)
{
    slong i;
    n_poly_t r;

    _n_poly_vec_mod_content(g, A, Alen, ctx);

    if (n_poly_is_one(g))
        return;

    n_poly_init(r);

    for (i = 0; i < Alen; i++)
        n_poly_mod_divrem(A + i, r, A + i, g, ctx);

    n_poly_clear(r);
}

/* TODO get rid of the need for this */
void n_polyun_mod_content(n_poly_t c, const n_polyun_t A, nmod_t ctx)
{
    slong i;

    n_poly_zero(c);
    for (i = 0; i < A->length; i++)
        n_poly_mod_gcd(c, c, A->coeffs + i, ctx);
}

