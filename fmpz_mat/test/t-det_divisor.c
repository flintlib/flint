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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    flint_rand_t state;
    long i;
    int result;

    printf("det_divisor....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        fmpz_t det, d, q, r;
        long m, bits;

        m = n_randint(state, 15);
        bits = 1 + n_randint(state, 50);

        fmpz_init(det);
        fmpz_init(d);
        fmpz_init(q);
        fmpz_init(r);
        fmpz_mat_init(A, m, m);

        if (i % 3 == 0 && m > 1)
        {
            /* Generate a nontrivial singular matrix */
            fmpz_mat_randrank(A, state, 1 + n_randint(state, m - 1), bits);
            fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));
        }
        else
        {
            fmpz_mat_randtest(A, state, bits);
        }

        fmpz_mat_det_divisor(d, A);
        fmpz_mat_det_bareiss(det, A);

        if (fmpz_is_zero(det) || fmpz_is_zero(d))
        {
            result = fmpz_equal(det, d);
        }
        else
        {
            fmpz_fdiv_qr(q, r, det, d);
            result = fmpz_is_zero(r) && (fmpz_sgn(d) > 0);
        }

        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mat_print_pretty(A), printf("\n");
            printf("det: ");  fmpz_print(det);    printf("\n");
            printf("d: "); fmpz_print(d); printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det);
        fmpz_clear(d);
        fmpz_clear(q);
        fmpz_clear(r);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
