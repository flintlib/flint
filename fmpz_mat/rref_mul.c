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

    Copyright (C) 2010-2012 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "perm.h"

/* checks matrix is in rref TODO are all these necessary here? */
int in_rref(const fmpz_mat_t A, const fmpz_t den, slong rank)
{
    slong i, j, k, prev_pivot;

    /* bottom should be zero */
    for (i = rank; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
                return 0;

    prev_pivot = -1;

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                /* pivot should have a higher column index than previous */
                if (j <= prev_pivot)
                    return 0;

                /* column should be 0 ... 0 1 0 ... 0 */
                for (k = 0; k < rank; k++)
                {
                    if (i == k && !fmpz_equal(fmpz_mat_entry(A, k, j), den))
                        return 0;
                    if (i != k && !fmpz_is_zero(fmpz_mat_entry(A, k, j)))
                        return 0;
                }

                prev_pivot = j;
                break;
            }
        }
    }

    return 1;
}

slong
fmpz_mat_rref_mul(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    flint_rand_t state;
    nmod_mat_t Amod;
    mp_limb_t p;
    slong i, j, m, n, rank, * pivs, * P;
    fmpz_mat_t B, C, D, E, E2, F, FD;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    flint_randinit(state);
    pivs = (slong *) flint_malloc(n * sizeof(slong));
    P = _perm_init(m);

    while (1)
    {
        p = n_randprime(state, NMOD_MAT_OPTIMAL_MODULUS_BITS, 1);
        nmod_mat_init(Amod, m, n, p);
        fmpz_mat_get_nmod_mat(Amod, A);

        rank = _nmod_mat_rref(Amod, pivs, P);

        nmod_mat_clear(Amod);

        fmpz_mat_init(B, rank, rank);
        fmpz_mat_init(C, rank, n - rank);

        /* set B to be the pivot columns and rows and C to be the non-pivot
           columns in the pivot rows */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
                fmpz_set(fmpz_mat_entry(B, i, j),
                        fmpz_mat_entry(A, P[i], pivs[j]));
            for (j = 0; j < n - rank; j++)
                fmpz_set(fmpz_mat_entry(C, i, j),
                        fmpz_mat_entry(A, P[i], pivs[rank + j]));
        }

        /* solve B*E2 = den*C */
        fmpz_mat_init(E2, rank, n - rank);
        if (n < 25)                 /* small matrices use solve */
        {
            if (!fmpz_mat_solve(E2, den, B, C))
            {
                flint_printf("Exception (fmpz_mat_rref_mul). "
                             "Singular input matrix for solve.");
                abort();
            }
        }
        else                        /* larger matrices use dixon */
        {
            fmpq_mat_t E2_q;
            if (!fmpz_mat_solve_dixon(E2, den, B, C))
            {
                flint_printf("Exception (fmpz_mat_rref_mul). "
                             "Singular input matrix for solve.");
                abort();
            }
            fmpq_mat_init(E2_q, rank, n - rank);
            fmpq_mat_set_fmpz_mat_mod_fmpz(E2_q, E2, den);
            fmpq_mat_get_fmpz_mat_matwise(E2, den, E2_q);
            fmpq_mat_clear(E2_q);
        }
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_init(E, rank, n);

        /* move columns of E2 and identity matrix into E so that it should be
           in rref */
        for (i = 0; i < rank; i++)
        {
            fmpz_set(fmpz_mat_entry(E, i, pivs[i]), den);
            for (j = 0; j < n - rank; j++)
                fmpz_set(fmpz_mat_entry(E, i, pivs[rank + j]),
                        fmpz_mat_entry(E2, i, j));
        }
        fmpz_mat_clear(E2);

        if (!in_rref(E, den, rank))
        {
            fmpz_mat_clear(E);
            continue;
        }

        /* set D to be the nullspace basis vector for E */
        fmpz_mat_init(D, n, n - rank);

        for (j = 0; j < n - rank; j++)
        {
            fmpz_set(fmpz_mat_entry(D, pivs[rank + j], j), den);
            for (i = 0; i < rank; i++)
                fmpz_neg(fmpz_mat_entry(D, pivs[i], j),
                        fmpz_mat_entry(E, i, pivs[rank + j]));
        }

        fmpz_mat_init(F, m - rank, n);

        for (i = 0; i < m - rank; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(F, i, j),
                        fmpz_mat_entry(A, P[rank + i], j));

        fmpz_mat_init(FD, m - rank, n - rank);
        fmpz_mat_mul(FD, F, D);
        fmpz_mat_clear(F);
        fmpz_mat_clear(D);

        /* if FD = 0 we have computed the right so stop, otherwise try a
           different p in the next iteration */
        if (fmpz_mat_is_zero(FD))
            break;

        fmpz_mat_clear(E);
        fmpz_mat_clear(FD);
    }

    /* write the entries of E  into R and zeroes at the bottom */
    for (i = 0; i < rank; i++)
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(R, i, j), fmpz_mat_entry(E, i, j));
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_zero(fmpz_mat_entry(R, i, j));

    fmpz_mat_clear(E);
    fmpz_mat_clear(FD);
    flint_randclear(state);
    flint_free(pivs);
    _perm_clear(P);

    return rank;
}

