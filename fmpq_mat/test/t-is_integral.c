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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("is_integral....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A;
        fmpz_mat_t B;
        long n, m;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);

        fmpz_mat_randtest(B, state, 10);
        fmpq_mat_set_fmpz_mat(A, B);

        if (!fmpq_mat_is_integral(A))
        {
            printf("FAIL\n");
            abort();
        }

        if (m && n)
        {
            long i, j;
            i = n_randint(state, m);
            j = n_randint(state, n);

            fmpq_set_si(fmpq_mat_entry(A, i, j), 1, 2);

            if (fmpq_mat_is_integral(A))
            {
                printf("FAIL\n");
                abort();
            }
        }


        fmpq_mat_clear(A);
        fmpz_mat_clear(B);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
