/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

#define E fmpz_poly_mat_entry

static inline void
fmpz_poly_addlow(fmpz_poly_t c, const fmpz_poly_t a,
    const fmpz_poly_t b, slong len)
{
    fmpz_poly_add(c, a, b);
    fmpz_poly_truncate(c, len);
}

void
fmpz_poly_mat_sqrlow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, slong len)
{
    slong n = A->r;

    if (n == 0)
        return;

    if (len < 1)
    {
        fmpz_poly_mat_zero(B);
        return;
    }

    if (n == 1)
    {
        fmpz_poly_sqrlow(E(B, 0, 0), E(A, 0, 0), len);
        return;
    }

    if (n == 2)
    {
        fmpz_poly_t t, u;

        fmpz_poly_init(t);
        fmpz_poly_init(u);

        fmpz_poly_addlow(t, E(A, 0, 0), E(A, 1, 1), len);
        fmpz_poly_mullow(u, E(A, 0, 1), E(A, 1, 0), len);

        fmpz_poly_sqrlow(E(B, 0, 0), E(A, 0, 0), len);
        fmpz_poly_addlow(E(B, 0, 0), E(B, 0, 0), u, len);

        fmpz_poly_sqrlow(E(B, 1, 1), E(A, 1, 1), len);
        fmpz_poly_addlow(E(B, 1, 1), E(B, 1, 1), u, len);

        fmpz_poly_mullow(E(B, 0, 1), E(A, 0, 1), t, len);
        fmpz_poly_mullow(E(B, 1, 0), E(A, 1, 0), t, len);

        fmpz_poly_clear(t);
        fmpz_poly_clear(u);
        return;
    }

    fmpz_poly_mat_mullow(B, A, A, len);
}
