/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "perm.h"

#define XX(ii,jj) fmpz_poly_mat_entry(X,(ii),(jj))
#define BB(ii,jj) fmpz_poly_mat_entry(B,(ii),(jj))
#define LU(ii,jj) fmpz_poly_mat_entry(FFLU,(ii),(jj))

void
fmpz_poly_mat_set_perm(fmpz_poly_mat_t X, const slong * perm,
    const fmpz_poly_mat_t B)
{
    if (X == B)
    {
        flint_throw(FLINT_ERROR, "(%s): Not implemented\n", __func__);
    }
    else
    {
        slong i, j;

        if (perm == NULL)
            flint_throw(FLINT_ERROR, "(%s): perm == NULL\n", __func__);

        for (i = 0; i < fmpz_poly_mat_nrows(B); i++)
            for (j = 0; j < fmpz_poly_mat_ncols(B); j++)
                fmpz_poly_set(fmpz_poly_mat_entry(X, i, j),
                              fmpz_poly_mat_entry(B, perm[i], j));
    }
}

void
fmpz_poly_mat_solve_fflu_precomp(fmpz_poly_mat_t X,
                    const slong * perm,
                    const fmpz_poly_mat_t FFLU, const fmpz_poly_mat_t B)
{
    fmpz_poly_t T;
    slong i, j, k, m, n;

    n = X->r;
    m = X->c;

    fmpz_poly_init(T);
    fmpz_poly_mat_set_perm(X, perm, B);

    for (k = 0; k < m; k++)
    {
        /* Fraction-free forward substitution */
        for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                fmpz_poly_mul(XX(j, k), XX(j, k), LU(i, i));
                fmpz_poly_mul(T, LU(j, i), XX(i, k));
                fmpz_poly_sub(XX(j, k), XX(j, k), T);
                if (i > 0)
                    fmpz_poly_div(XX(j, k), XX(j, k), LU(i-1, i-1));
            }
        }

        /* Fraction-free back substitution */
        for (i = n - 2; i >= 0; i--)
        {
            fmpz_poly_mul(XX(i, k), XX(i, k), LU(n-1, n-1));
            for (j = i + 1; j < n; j++)
            {
                fmpz_poly_mul(T, XX(j, k), LU(i, j));
                fmpz_poly_sub(XX(i, k), XX(i, k), T);
            }
            fmpz_poly_div(XX(i, k), XX(i, k), LU(i, i));
        }
    }

    fmpz_poly_clear(T);
}
