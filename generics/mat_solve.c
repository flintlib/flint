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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"
#include "perm.h"


void
elem_mat_set_perm(elem_mat_t X, const long * perm, const elem_mat_t B, const ring_t ring)
{
    long c, i, j, m, n;
    const ring_struct * ering = RING_PARENT(ring);

    m = elem_mat_nrows(B, ring);
    n = elem_mat_ncols(B, ring);

    if (X == B)
    {
        elem_ptr tmp = _elem_vec_init(n, ering);

        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
                elem_set(INDEX(tmp, i, ering->size), elem_mat_entry(B, perm[i], c, ring), ering);
            for (i = 0; i < n; i++)
                elem_set(elem_mat_entry(X, i, c, ring), INDEX(tmp, i, ering->size), ering);
        }

        _elem_vec_clear(tmp, n, ering);
    }
    else
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                elem_set(elem_mat_entry(X, i, j, ring), elem_mat_entry(B, perm[i], j, ring), ering);
    }
}

#define XX(ii,jj) elem_mat_entry(X,(ii),(jj),ring)
#define BB(ii,jj) elem_mat_entry(B,(ii),(jj),ring)
#define LU(ii,jj) elem_mat_entry(FFLU,(ii),(jj),ring)

void
elem_mat_solve_fflu_precomp(elem_mat_t X,
                    const long * perm,
                    const elem_mat_t FFLU, const elem_mat_t B, const ring_t ring)
{
    elem_ptr T;
    long i, j, k, m, n;
    const ring_struct * ering = RING_PARENT(ring);

    n = X->r;
    m = X->c;

    ELEM_TMP_INIT(T, ering);
    elem_mat_set_perm(X, perm, B, ring);

    for (k = 0; k < m; k++)
    {
        /* Fraction-free forward substitution */
        for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                elem_mul(XX(j, k), XX(j, k), LU(i, i), ering);
                elem_mul(T, LU(j, i), XX(i, k), ering);
                elem_sub(XX(j, k), XX(j, k), T, ering);
                if (i > 0)
                {
                    elem_divexact(XX(j, k), XX(j, k), LU(i-1, i-1), ering);
                }
            }
        }

        /* Fraction-free back substitution */
        for (i = n - 2; i >= 0; i--)
        {
            elem_mul(XX(i, k), XX(i, k), LU(n-1, n-1), ering);
            for (j = i + 1; j < n; j++)
            {
                elem_mul(T, XX(j, k), LU(i, j), ering);
                elem_sub(XX(i, k), XX(i, k), T, ering);
            }
            elem_divexact(XX(i, k), XX(i, k), LU(i, i), ering);
        }
    }

    ELEM_TMP_CLEAR(T, ering);
}

int
elem_mat_solve(elem_mat_t X, elem_ptr den, const elem_mat_t A, const elem_mat_t B, const ring_t ring)
{
    elem_mat_t LU;
    long dim, *perm;
    int result;

    if (elem_mat_is_empty(B, ring))
    {
        elem_one(den, RING_PARENT(ring));
        return 1;
    }

    dim = elem_mat_nrows(A, ring);
    perm = _perm_init(dim);
    elem_mat_init(LU, dim, dim, ring);
    elem_mat_set(LU, A, ring);

    result = (elem_mat_fflu(LU, den, perm, LU, 1, ring) == dim);

    if (result)
        elem_mat_solve_fflu_precomp(X, perm, LU, B, ring);
    else
        elem_zero(den, RING_PARENT(ring));

    _perm_clear(perm);
    elem_mat_clear(LU, ring);
    return result;
}

