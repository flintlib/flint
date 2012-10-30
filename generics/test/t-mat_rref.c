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

#include <stdlib.h>
#include <stdio.h>
#include "generics.h"
#include "perm.h"


/* checks that the rref has the right form */
int check_rref(const elem_mat_t A, elem_srcptr den, long rank, const ring_t ring)
{
    long i, j, k, prev_pivot;
    ring_struct * ering = RING_PARENT(ring);

    /* bottom should be zero */
    for (i = rank; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (!elem_is_zero(elem_mat_entry(A, i, j, ring), ering))
                return 0;

    prev_pivot = -1;

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (!elem_is_zero(elem_mat_entry(A, i, j, ring), ering))
            {
                /* pivot should have a higher column index than previous */
                if (j <= prev_pivot)
                    return 0;

                /* column should be 0 ... 0 1 0 ... 0 */
                for (k = 0; k < rank; k++)
                {
                    if (i == k && !elem_equal(elem_mat_entry(A, k, j, ring), den, ering))
                        return 0;
                    if (i != k && !elem_is_zero(elem_mat_entry(A, k, j, ring), ering))
                        return 0;
                }

                prev_pivot = j;
                break;
            }
        }
    }

    return 1;
}


int
main(void)
{
    long iter;
    flint_rand_t state;

    printf("mat_rref....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        elem_mat_t A, R, B, R2;
        elem_poly_t den, c, den2;
        ring_t ZZ, ZZx, MM;
        long j, k, m, n, rank1, rank2;
        long size[3] = {5, 5, 5};
        long *perm;
        int equal;

        m = n_randint(state, 7);
        n = n_randint(state, 7);

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_mat(MM, ZZx);

        elem_mat_init(A, m, n, MM);
        elem_mat_init(R, m, n, MM);
        elem_mat_init(B, 2 * m, n, MM);
        elem_mat_init(R2, 2 * m, n, MM);

        elem_init(c, ZZx);
        elem_init(den, ZZx);
        elem_init(den2, ZZx);

        perm = _perm_init(2 * m);

        elem_mat_randtest(A, state, size, MM);

        rank1 = elem_mat_rref(R, den, A, MM);

        check_rref(R, den, rank1, MM);

        /* Concatenate the original matrix with the rref,
            scramble the rows, and check that the rref is the same */
        _perm_randtest(perm, 2 * m, state);

        for (j = 0; j < m; j++)
        {
            elem_randtest_not_zero(c, state, size, ZZx);
            for (k = 0; k < n; k++)
                elem_mul(elem_mat_entry(B, perm[j], k, MM),
                    elem_mat_entry(A, j, k, MM), c, ZZx);
        }

        for (j = 0; j < m; j++)
        {
            elem_randtest_not_zero(c, state, size, ZZx);
            for (k = 0; k < n; k++)
                elem_mul(elem_mat_entry(B, perm[m + j], k, MM),
                    elem_mat_entry(R, j, k, MM), c, ZZx);
        }

        rank2 = elem_mat_rref(R2, den2, B, MM);
        equal = (rank1 == rank2);

        if (equal)
        {
            elem_mat_scalar_mul(R, R, den2, MM);
            elem_mat_scalar_mul(R2, R2, den, MM);

            for (j = 0; j < rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal &&
                        elem_equal(elem_mat_entry(R, j, k, MM), elem_mat_entry(R2, j, k, MM), ZZx);
            for (j = rank2; j < 2 * rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal && elem_is_zero(elem_mat_entry(R2, j, k, MM), ZZx);
        }

        if (!equal)
        {
            printf("FAIL (rank1 = %ld, rank2 = %ld)!\n", rank1, rank2);
            elem_print(A, MM); printf("\n\n");
            elem_print(R, MM); printf("\n\n");
            elem_print(R2, MM); printf("\n\n");
            abort();
        }

        elem_clear(c, ZZx);
        elem_clear(den, ZZx);
        elem_clear(den2, ZZx);

        _perm_clear(perm);

        elem_mat_clear(A, MM);
        elem_mat_clear(R, MM);
        elem_mat_clear(B, MM);
        elem_mat_clear(R2, MM);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

