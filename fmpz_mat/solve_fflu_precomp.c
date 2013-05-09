/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

#define XX(ii,jj) fmpz_mat_entry(X,(ii),(jj))
#define BB(ii,jj) fmpz_mat_entry(B,(ii),(jj))
#define LU(ii,jj) fmpz_mat_entry(FFLU,(ii),(jj))

void
fmpz_mat_set_perm(fmpz_mat_t X, const len_t * perm, const fmpz_mat_t B)
{
    if (X == B)
    {
        /* Not implemented */
        abort();
    }
    else
    {
        len_t i, j;

        if (perm == NULL)
            abort();

        for (i = 0; i < fmpz_mat_nrows(B); i++)
            for (j = 0; j < fmpz_mat_ncols(B); j++)
                fmpz_set(fmpz_mat_entry(X, i, j),
                         fmpz_mat_entry(B, perm[i], j));
    }
}

void
fmpz_mat_solve_fflu_precomp(fmpz_mat_t X,
                    const len_t * perm,
                    const fmpz_mat_t FFLU, const fmpz_mat_t B)
{
    fmpz_t T;
    len_t i, j, k, m, n;

    n = X->r;
    m = X->c;

    fmpz_init(T);
    fmpz_mat_set_perm(X, perm, B);

    for (k = 0; k < m; k++)
    {
        /* Fraction-free forward substitution */
        for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                fmpz_mul(XX(j, k), XX(j, k), LU(i, i));
                fmpz_mul(T, LU(j, i), XX(i, k));
                fmpz_sub(XX(j, k), XX(j, k), T);
                if (i > 0)
                    fmpz_divexact(XX(j, k), XX(j, k), LU(i-1, i-1));
            }
        }

        /* Fraction-free back substitution */
        for (i = n - 2; i >= 0; i--)
        {
            fmpz_mul(XX(i, k), XX(i, k), LU(n-1, n-1));
            for (j = i + 1; j < n; j++)
            {
                fmpz_mul(T, XX(j, k), LU(i, j));
                fmpz_sub(XX(i, k), XX(i, k), T);
            }
            fmpz_divexact(XX(i, k), XX(i, k), LU(i, i));
        }
    }

    fmpz_clear(T);
}
