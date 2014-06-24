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

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"

void fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong j, j2, i, k, l;
    fmpz_t u, v, d, r2d, r1d, q, b;

    fmpz_mat_set(H, A);
    /* old preprocessing step TODO remove when completely sure nothing is broken
    fmpz_init(det);
    fmpz_mat_init(M, n, n);
    fmpz_mat_one(M);
    ensure all principal minors are non-singular
    for (i = 0; i < n; i++)
    {
        for (i2 = 0; i2 < i; i2++)
            fmpz_set(fmpz_mat_entry(M, i2, i), fmpz_mat_entry(H, i2, i));
        fmpz_zero(det);
        for (i2 = i; i2 < n; i2++)
        {
            for (j = 0; j <= i; j++)
                fmpz_set(fmpz_mat_entry(M, i, j), fmpz_mat_entry(H, i2, j));
            fmpz_mat_det(det, M);
            if (!fmpz_is_zero(det))
                break;
        }
        if (i2 == n)
        {
            flint_printf("Exception (fmpz_mat_hnf_minors). Input matrix was singular.\n");
            abort();
        }
        fmpz_mat_swap_rows(H, NULL, i2, i);
    }
    fmpz_clear(det);
    fmpz_mat_clear(M); */
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(d);
    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(q);
    fmpz_init(b);
    /* put the kth principal minor in HNF */
    for (k = 0, l = A->r - 1; k < A->c; k++)
    {
        for (j = 0; j < k; j++)
        {
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, j, j), fmpz_mat_entry(H, k, j));
            fmpz_divexact(r1d, fmpz_mat_entry(H, j, j), d);
            fmpz_divexact(r2d, fmpz_mat_entry(H, k, j), d);
            for (j2 = j; j2 < A->c; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, j, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, k, j2));
                fmpz_mul(fmpz_mat_entry(H, k, j2), r1d, fmpz_mat_entry(H, k, j2));
                fmpz_submul(fmpz_mat_entry(H, k, j2), r2d, fmpz_mat_entry(H, j, j2));
                fmpz_set(fmpz_mat_entry(H, j, j2), b);
            }
            /* ensure H_j,j is positive
            if (fmpz_sgn(fmpz_mat_entry(H, j, j)) < 0)
            {
                for (j2 = j; j2 < A->c; j2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, j, j2),
                            fmpz_mat_entry(H, j, j2));
                }
            } */
        }
        /* if H_k,k is zero we swap row k for some other row (starting with the last) */
        if (fmpz_is_zero(fmpz_mat_entry(H, k, k)))
        {
            fmpz_mat_swap_rows(H, NULL, k, l);
            l--;
            k--;
            continue;
        }
        /* ensure H_k,k is positive */
        if (fmpz_sgn(fmpz_mat_entry(H, k, k)) < 0)
        {
            for (j = k; j < A->c; j++)
            {
                fmpz_neg(fmpz_mat_entry(H, k, j),
                        fmpz_mat_entry(H, k, j));
            }
        }
        /* reduce above diagonal elements of each row i */
        for (i = k - 1; i >= 0; i--)
        {
            for (j = i + 1; j <= k; j++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                        fmpz_mat_entry(H, j, j));
                for (j2 = j; j2 < A->c; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                            fmpz_mat_entry(H, j, j2));
                }
            }
        }
        l = A->r - 1;
    }
    /* reduce final rows */
    for (k = A->c; k < A->r; k++)
    {
        for (j = 0; j < A->c; j++)
        {
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, j, j), fmpz_mat_entry(H, k, j));
            fmpz_divexact(r1d, fmpz_mat_entry(H, j, j), d);
            fmpz_divexact(r2d, fmpz_mat_entry(H, k, j), d);
            for (j2 = j; j2 < A->c; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, j, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, k, j2));
                fmpz_mul(fmpz_mat_entry(H, k, j2), r1d, fmpz_mat_entry(H, k, j2));
                fmpz_submul(fmpz_mat_entry(H, k, j2), r2d, fmpz_mat_entry(H, j, j2));
                fmpz_set(fmpz_mat_entry(H, j, j2), b);
            }
        }
        /* reduce above diagonal elements of each row i */
        for (i = A->c - 1; i >= 0; i--)
        {
            for (j = i + 1; j < A->c; j++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                        fmpz_mat_entry(H, j, j));
                for (j2 = j; j2 < A->c; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                            fmpz_mat_entry(H, j, j2));
                }
            }
        }
    }
    fmpz_clear(b);
    fmpz_clear(q);
    fmpz_clear(r2d);
    fmpz_clear(r1d);
    fmpz_clear(d);
    fmpz_clear(v);
    fmpz_clear(u);
}
