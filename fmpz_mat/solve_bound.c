/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_solve_bound(fmpz_t N, fmpz_t D,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong i, j, m, n;
    fmpz_t t, u;

    m = B->r;
    n = B->c;

    fmpz_mat_det_bound(D, A);

    fmpz_init(t);
    fmpz_init(u);

    fmpz_zero(t);

    /* Largest column norm of B */
    for (j = 0; j < n; j++)
    {
        fmpz_zero(u);
        for (i = 0; i < m; i++)
            fmpz_addmul(u, fmpz_mat_entry(B, i, j), fmpz_mat_entry(B, i, j));
        if (fmpz_cmp(t, u) < 0)
            fmpz_set(t, u);
    }

    fmpz_sqrtrem(t, u, t);
    if (!fmpz_is_zero(u))
        fmpz_add_ui(t, t, UWORD(1));

    fmpz_mul(N, D, t);

    fmpz_clear(t);
    fmpz_clear(u);
}
