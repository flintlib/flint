/*
    Copyright (C) 2011, 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fmpz_mat.h"
#include "nmod_poly_mat.h"

#define KS_MIN_DIM 10
#define INTERPOLATE_MIN_DIM 80
#define KS_MAX_LENGTH 128

void
nmod_poly_mat_sqr(nmod_poly_mat_t C, const nmod_poly_mat_t A)
{
    ulong dim = A->r;

    if (dim < KS_MIN_DIM)
    {
        nmod_poly_mat_sqr_classical(C, A);
    }
    else
    {
        ulong Alen;
        ulong mod = nmod_poly_mat_modulus(A);

        Alen = nmod_poly_mat_max_length(A);

        if ((FLINT_BIT_COUNT(mod) > FLINT_BITS / 4)
            && (dim > INTERPOLATE_MIN_DIM + n_sqrt(Alen))
            && (mod >= 2 * Alen - 1) && n_is_prime(mod))
            nmod_poly_mat_sqr_interpolate(C, A);

        if (Alen > KS_MAX_LENGTH)
            nmod_poly_mat_sqr_classical(C, A);
        else
            nmod_poly_mat_sqr_KS(C, A);
    }
}

static inline void
nmod_poly_sqr(nmod_poly_t y, const nmod_poly_t x)
{
    nmod_poly_mul(y, x, x);
}

#define E nmod_poly_mat_entry

void
nmod_poly_mat_sqr_classical(nmod_poly_mat_t B, const nmod_poly_mat_t A)
{
    slong n = A->r;

    if (n == 0)
        return;

    if (n == 1)
    {
        nmod_poly_sqr(E(B, 0, 0), E(A, 0, 0));
        return;
    }

    if (n == 2)
    {
        nmod_poly_t t, u;

        nmod_poly_init(t, nmod_poly_mat_modulus(A));
        nmod_poly_init(u, nmod_poly_mat_modulus(A));

        nmod_poly_add(t, E(A, 0, 0), E(A, 1, 1));
        nmod_poly_mul(u, E(A, 0, 1), E(A, 1, 0));

        nmod_poly_sqr(E(B, 0, 0), E(A, 0, 0));
        nmod_poly_add(E(B, 0, 0), E(B, 0, 0), u);

        nmod_poly_sqr(E(B, 1, 1), E(A, 1, 1));
        nmod_poly_add(E(B, 1, 1), E(B, 1, 1), u);

        nmod_poly_mul(E(B, 0, 1), E(A, 0, 1), t);
        nmod_poly_mul(E(B, 1, 0), E(A, 1, 0), t);

        nmod_poly_clear(t);
        nmod_poly_clear(u);
        return;
    }

    nmod_poly_mat_mul_classical(B, A, A);
}

void
nmod_poly_mat_sqr_interpolate(nmod_poly_mat_t C, const nmod_poly_mat_t A)
{
    nmod_poly_mat_mul_interpolate(C, A, A);
}

void
nmod_poly_mat_sqr_KS(nmod_poly_mat_t B, const nmod_poly_mat_t A)
{
    nmod_poly_mat_mul_KS(B, A, A);
}
