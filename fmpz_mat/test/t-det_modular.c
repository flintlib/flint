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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    fmpz_mat_t A;
    flint_rand_t state;
    long i, m;

    fmpz_t det1, det2;

    printf("det_modular....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int proved = n_randlimb(state) % 2;
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        fmpz_init(det1);
        fmpz_init(det2);

        fmpz_mat_randtest(A, state, 1+n_randint(state,200));

        fmpz_mat_det_bareiss(det1, A);
        fmpz_mat_det_modular(det2, A, proved);

        if (!fmpz_equal(det1, det2))
        {
            printf("FAIL:\n");
            printf("different determinants!\n");
            fmpz_mat_print_pretty(A), printf("\n");
            printf("det1: "), fmpz_print(det1), printf("\n");
            printf("det2: "), fmpz_print(det2), printf("\n");
            abort();
        }

        fmpz_clear(det1);
        fmpz_clear(det2);
        fmpz_mat_clear(A);
    }

    for (i = 0; i < 10000; i++)
    {
        int proved = n_randlimb(state) % 2;
        m = 2 + n_randint(state, 10);
        fmpz_mat_init(A, m, m);
        fmpz_init(det2);

        fmpz_mat_randrank(A, state, 1+n_randint(state, m - 1),
                                    1+n_randint(state, 10));
        fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));

        fmpz_mat_det_modular(det2, A, proved);
        if (!fmpz_is_zero(det2))
        {
            printf("FAIL:\n");
            printf("expected zero determinant!\n");
            fmpz_mat_print_pretty(A), printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
