/*
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_sqr_classical(arb_mat_t B, const arb_mat_t A, slong prec)
{
    slong n, i, j, k;
    arb_t p, s;

    n = arb_mat_nrows(A);

    if (arb_mat_ncols(A) != n || arb_mat_nrows(B) != n || arb_mat_ncols(B) != n)
    {
        flint_throw(FLINT_ERROR, "arb_mat_sqr: incompatible dimensions\n");
    }

    if (n == 0)
        return;

    if (n == 1)
    {
        arb_mul(arb_mat_entry(B, 0, 0),
                  arb_mat_entry(A, 0, 0),
                  arb_mat_entry(A, 0, 0), prec);
        return;
    }

    if (A == B)
    {
        arb_mat_t T;
        arb_mat_init(T, n, n);
        arb_mat_sqr_classical(T, A, prec);
        arb_mat_swap(T, B);
        arb_mat_clear(T);
        return;
    }

    arb_init(p);
    arb_init(s);

    /* contribution of diagonal of A to diagonal of B */
    for (i = 0; i < n; i++)
    {
        arb_mul(arb_mat_entry(B, i, i),
                  arb_mat_entry(A, i, i),
                  arb_mat_entry(A, i, i), prec);
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            /* contribution of off-diagonal of A to diagonal of B */
            arb_mul(p, arb_mat_entry(A, i, j), arb_mat_entry(A, j, i), prec);
            arb_add(arb_mat_entry(B, i, i), arb_mat_entry(B, i, i), p, prec);
            arb_add(arb_mat_entry(B, j, j), arb_mat_entry(B, j, j), p, prec);

            /* contribution of diagonal of A to off-diagonal of B */
            arb_add(s, arb_mat_entry(A, i, i), arb_mat_entry(A, j, j), prec);
            arb_mul(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), s, prec);
            arb_mul(arb_mat_entry(B, j, i), arb_mat_entry(A, j, i), s, prec);
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
                        arb_addmul(arb_mat_entry(B, i, j),
                                     arb_mat_entry(A, i, k),
                                     arb_mat_entry(A, k, j), prec);
                    }
                }
            }
        }
    }

    arb_clear(p);
    arb_clear(s);
}

