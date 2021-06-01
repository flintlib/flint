/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


slong _n_poly_vec_max_degree(const n_poly_struct * A, slong Alen)
{
    slong i, len = 0;
    for (i = 0; i < Alen; i++)
        len = FLINT_MAX(len, A[i].length);
    return len - 1;
}


void _n_poly_vec_mul_nmod_intertible(
    n_poly_struct * A,
    slong Alen,
    mp_limb_t c,
    nmod_t ctx)
{
    slong i;

    FLINT_ASSERT(n_gcd(c, ctx.n) == 1);

    if (c == 1)
        return;

    for (i = 0; i < Alen; i++)
        _n_poly_mod_scalar_mul_nmod_inplace(A + i, c, ctx);
}


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

void _n_poly_vec_mod_divexact_poly(
    n_poly_struct * A,
    slong Alen,
    const n_poly_t g,
    nmod_t ctx)
{
    slong i;
    n_poly_t r;

    if (n_poly_is_one(g))
        return;

    n_poly_init(r);

    for (i = 0; i < Alen; i++)
    {
        n_poly_mod_divrem(A + i, r, A + i, g, ctx);
        FLINT_ASSERT(n_poly_is_zero(r));
    }

    n_poly_clear(r);
}


void _n_poly_vec_mod_remove_content(
    n_poly_t g,
    n_poly_struct * A,
    slong Alen,
    nmod_t ctx)
{
    _n_poly_vec_mod_content(g, A, Alen, ctx);
    _n_poly_vec_mod_divexact_poly(A, Alen, g, ctx);
}

