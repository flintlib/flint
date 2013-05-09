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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

int
main(void)
{
    flint_rand_t state;
    int iter;

    printf("one/is_one....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_mat_t A;
        len_t m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_randtest(A, state, n_randint(state, 5),
            n_randint(state, 100));
        fmpz_poly_mat_one(A);

        if (!fmpz_poly_mat_is_one(A))
        {
            printf("FAIL: expected matrix to be one\n");
            abort();
        }

        if (m > 0 && n > 0)
        {
            m = n_randint(state, m);
            n = n_randint(state, n);

            if (m != n)
                fmpz_poly_randtest_not_zero(fmpz_poly_mat_entry(A, m, n),
                    state, 5, 5);
            else
                do { fmpz_poly_randtest_not_zero(fmpz_poly_mat_entry(A, m, n),
                    state, 5, 5); }
                while (fmpz_poly_is_one(fmpz_poly_mat_entry(A, m, n)));

            if (fmpz_poly_mat_is_one(A))
            {
                printf("FAIL: expected matrix not to be one\n");
                abort();
            }
        }

        fmpz_poly_mat_clear(A);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
