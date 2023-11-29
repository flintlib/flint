/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

static void
_arb_sqr(arb_t dest, const arb_t src, slong prec)
{
    arb_mul(dest, src, src, prec);
}

void
arb_mat_inv_cho_precomp(arb_mat_t X, const arb_mat_t L, slong prec)
{
    slong n;

    if (arb_mat_nrows(X) != arb_mat_nrows(L) ||
        arb_mat_ncols(X) != arb_mat_ncols(L))
    {
        flint_throw(FLINT_ERROR, "arb_mat_inv_cho_precomp: incompatible dimensions\n");
    }

    if (arb_mat_is_empty(L))
        return;

    n = arb_mat_nrows(L);

    if (n == 1)
    {
        arb_inv(arb_mat_entry(X, 0, 0), arb_mat_entry(L, 0, 0), prec);
        _arb_sqr(arb_mat_entry(X, 0, 0), arb_mat_entry(X, 0, 0), prec);
        return;
    }

    if (X == L)
    {
        flint_throw(FLINT_ERROR, "arb_mat_inv_cho_precomp: unsupported aliasing\n");
    }

    /* invert a 2x2 or larger matrix given its L * L^T decomposition */
    {
        slong i, j, k;
        arb_struct *s;
        arb_mat_zero(X);
        s = _arb_vec_init(n);
        for (i = 0; i < n; i++)
        {
            arb_inv(s + i, arb_mat_entry(L, i, i), prec);
        }
        for (j = n-1; j >= 0; j--)
        {
            for (i = j; i >= 0; i--)
            {
                if (i == j)
                {
                    arb_set(arb_mat_entry(X, i, j), s + i);
                }
                else
                {
                    arb_zero(arb_mat_entry(X, i, j));
                }
                for (k = i + 1; k < n; k++)
                {
                    arb_submul(arb_mat_entry(X, i, j),
                               arb_mat_entry(L, k, i),
                               arb_mat_entry(X, k, j), prec);
                }
                arb_div(arb_mat_entry(X, i, j),
                        arb_mat_entry(X, i, j),
                        arb_mat_entry(L, i, i), prec);
                arb_set(arb_mat_entry(X, j, i),
                        arb_mat_entry(X, i, j));
            }
        }
        _arb_vec_clear(s, n);
    }
}
