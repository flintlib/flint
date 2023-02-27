/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod)
{
    slong i;

    if (A->length < 1)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_mod_is_canonical(A->coeffs + i, mod))
            return 0;
        if (i + 1 == A->length && n_poly_is_zero(A->coeffs + i))
            return 0;
    }

    return 1;
}


void n_bpoly_scalar_mul_nmod(n_bpoly_t A, mp_limb_t c, nmod_t ctx)
{
    slong i;

    if (c <= 1)
    {
        if (c == 0)
            A->length = 0;
        return;
    }

    for (i = 0; i < A->length; i++)
        _n_poly_mod_scalar_mul_nmod_inplace(A->coeffs + i, c, ctx);
}


void n_bpoly_mod_content_last(n_poly_t g, const n_bpoly_t A, nmod_t ctx)
{
    slong i;

    n_poly_zero(g);
    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_gcd(g, g, A->coeffs + i, ctx);
        if (n_poly_degree(g) == 0)
            break;
    }
}

void n_bpoly_mod_divexact_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx)
{
    slong i;
    n_poly_struct * t;

    if (b->length == 1)
    {
        if (b->coeffs[0] != 1)
            n_bpoly_scalar_mul_nmod(A, nmod_inv(b->coeffs[0], ctx), ctx);
        return;
    }

    n_bpoly_fit_length(A, A->length + 1);

    t = A->coeffs + A->length;

    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i].length < 1)
            continue;

        n_poly_mod_div(t, A->coeffs + i, b, ctx);
        n_poly_swap(A->coeffs + i, t);
    }
}


void n_bpoly_mod_mul_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx)
{
    slong i;
    n_poly_struct * t;

    if (n_poly_is_one(b))
        return;

    n_bpoly_fit_length(A, A->length + 1);

    t = A->coeffs + A->length;

    for (i = 0; i < A->length; i++)
    {
        if (n_poly_is_zero(A->coeffs + i))
            continue;

        n_poly_mod_mul(t, A->coeffs + i, b, ctx);
        n_poly_swap(A->coeffs + i, t);
    }
}

void n_bpoly_mod_derivative_gen0(
    n_bpoly_t A,
    const n_bpoly_t B,
    nmod_t ctx)
{
    slong i;

    FLINT_ASSERT(A != B);

    if (B->length < 2)
    {
        n_bpoly_zero(A);
        return;
    }

    n_bpoly_fit_length(A, B->length - 1);

    for (i = 1; i < B->length; i++)
        n_poly_mod_scalar_mul_ui(A->coeffs + i - 1, B->coeffs + i, i, ctx);

    A->length = B->length - 1;
    n_bpoly_normalise(A);
}
