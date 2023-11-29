/*
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_sqr_classical(acb_mat_t B, const acb_mat_t A, slong prec)
{
    slong n, i, j, k;
    acb_t p, s;

    n = acb_mat_nrows(A);

    if (acb_mat_ncols(A) != n || acb_mat_nrows(B) != n || acb_mat_ncols(B) != n)
    {
        flint_throw(FLINT_ERROR, "acb_mat_sqr: incompatible dimensions\n");
    }

    if (n == 0)
        return;

    if (n == 1)
    {
        acb_mul(acb_mat_entry(B, 0, 0),
                  acb_mat_entry(A, 0, 0),
                  acb_mat_entry(A, 0, 0), prec);
        return;
    }

    if (A == B)
    {
        acb_mat_t T;
        acb_mat_init(T, n, n);
        acb_mat_sqr_classical(T, A, prec);
        acb_mat_swap(T, B);
        acb_mat_clear(T);
        return;
    }

    acb_init(p);
    acb_init(s);

    /* contribution of diagonal of A to diagonal of B */
    for (i = 0; i < n; i++)
    {
        acb_mul(acb_mat_entry(B, i, i),
                  acb_mat_entry(A, i, i),
                  acb_mat_entry(A, i, i), prec);
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            /* contribution of off-diagonal of A to diagonal of B */
            acb_mul(p, acb_mat_entry(A, i, j), acb_mat_entry(A, j, i), prec);
            acb_add(acb_mat_entry(B, i, i), acb_mat_entry(B, i, i), p, prec);
            acb_add(acb_mat_entry(B, j, j), acb_mat_entry(B, j, j), p, prec);

            /* contribution of diagonal of A to off-diagonal of B */
            acb_add(s, acb_mat_entry(A, i, i), acb_mat_entry(A, j, j), prec);
            acb_mul(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), s, prec);
            acb_mul(acb_mat_entry(B, j, i), acb_mat_entry(A, j, i), s, prec);
        }
    }

    /* contribution of off-diagonal of A to off-diagonal of B */
    if (n > 2)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                for (k = 0; k < n; k++)
                {
                    if (i != j && j != k && k != i)
                    {
                        acb_addmul(acb_mat_entry(B, i, j),
                                     acb_mat_entry(A, i, k),
                                     acb_mat_entry(A, k, j), prec);
                    }
                }
            }
        }
    }

    acb_clear(p);
    acb_clear(s);
}

