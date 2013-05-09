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
    len_t i, m;

    fmpz_t det, bound;

    printf("det_bound....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        fmpz_init(det);
        fmpz_init(bound);

        fmpz_mat_randtest(A, state, 1+n_randint(state,200));

        fmpz_mat_det(det, A);
        fmpz_mat_det_bound(bound, A);

        if (fmpz_cmp(det, bound) > 0)
        {
            printf("FAIL:\n");
            printf("bound too small!\n");
            fmpz_mat_print_pretty(A), printf("\n");
            printf("det: "), fmpz_print(det), printf("\n");
            printf("bound: "), fmpz_print(bound), printf("\n");
            abort();
        }

        fmpz_clear(det);
        fmpz_clear(bound);
        fmpz_mat_clear(A);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
