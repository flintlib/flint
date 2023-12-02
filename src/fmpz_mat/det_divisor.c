/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq.h"

void
fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A)
{
    fmpz_mat_t X, B;
    fmpz_t t, u, v, mod;
    slong i, n;
    int success;

    n = A->r;

    fmpz_mat_init(B, n, 1);
    fmpz_mat_init(X, n, 1);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(mod);

    /* Create a "random" vector */
    for (i = 0; i < n; i++)
    {
        fmpz_set_si(fmpz_mat_entry(B, i, 0), 2*(i % 2) - 1);
    }

    success = fmpz_mat_solve_dixon(X, mod, A, B);

    if (success)
    {
        fmpz_one(d);
        for (i = 0; i < n; i++)
        {
            fmpz_mul(t, d, fmpz_mat_entry(X, i, 0));
            fmpz_fdiv_qr(u, t, t, mod);
            if (!_fmpq_reconstruct_fmpz(u, v, t, mod))
            {
                flint_throw(FLINT_ERROR, "(fmpz_mat_det_divisor): Rational reconstruction failed.\n");
            }

            fmpz_mul(d, v, d);
        }
    }
    else
    {
        fmpz_zero(d);
    }

    fmpz_mat_clear(B);
    fmpz_mat_clear(X);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(mod);
}
