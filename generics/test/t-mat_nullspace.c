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

int
main(void)
{
    long iter;
    flint_rand_t state;

    printf("mat_nullspace....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        elem_mat_t A, N, AN;
        ring_t ZZ, ZZx, MM;
        long m, n, rank, nullity;
        long size[3] = {5, 5, 5};

        m = n_randint(state, 7);
        n = n_randint(state, 7);

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_mat(MM, ZZx);

        elem_mat_init(A, m, n, MM);
        elem_mat_init(N, n, n, MM);
        elem_mat_init(AN, m, n, MM);

        elem_mat_randtest(A, state, size, MM);

        rank = elem_mat_rank(A, MM);
        nullity = elem_mat_nullspace(N, A, MM);

        if (nullity + rank != n)
        {
            printf("FAIL: wrong nullity!\n");
            printf("rank = %ld\n", rank);
            printf("nullity = %ld\n", nullity);
            elem_mat_print(A, MM);
            printf("\n");
            elem_mat_print(N, MM);
            printf("\n");
            abort();
        }

        if (elem_mat_rank(N, MM) != nullity)
        {
            printf("FAIL: wrong rank(N) != nullity!\n");
            abort();
        }

        elem_mat_mul(AN, A, N, MM);

        if (!elem_mat_is_zero(AN, MM))
        {
            printf("FAIL: A * N != 0\n");
            abort();
        }

        elem_mat_clear(A, MM);
        elem_mat_clear(N, MM);
        elem_mat_clear(AN, MM);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

