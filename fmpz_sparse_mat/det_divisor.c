/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "fmpq.h"

void
fmpz_sparse_mat_det_divisor(fmpz_t d, const fmpz_sparse_mat_t M)
{
    int success;
    slong i;
    fmpz_mat_t X, B;
    fmpz_t t, u, v, mod;
    if (M->r != M->c) {fmpz_zero(d); return;}

    fmpz_mat_init(X, M->c, 1);
    fmpz_mat_init(B, M->r, 1);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(mod);

    /* Create a "random" vector */
    for (i = 0; i < M->r; i++)
        fmpz_set_si(fmpz_mat_entry(B, i, 0), 2*(i % 2) - 1);

    success = fmpz_sparse_mat_solve_dixon(X, mod, M, B);
    if (success)
    {
        fmpz_one(d);
        for (i = 0; i < M->r; i++)
        {
            fmpz_mul(t, d, fmpz_mat_entry(X, i, 0));
            fmpz_fdiv_qr(u, t, t, mod);
            if (!_fmpq_reconstruct_fmpz(u, v, t, mod))
            {
                flint_printf("Exception (fmpz_mat_det_divisor): "
                       "Rational reconstruction failed.\n");
                flint_abort();
            }

            fmpz_mul(d, v, d);
        }
    }
    else fmpz_zero(d);

    fmpz_mat_clear(X);
    fmpz_mat_clear(B);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(mod);
}
