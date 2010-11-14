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
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    nmod_mat_t A, B, C, I;
    long i, j, m, r;
    mp_limb_t mod;
    int result;

    printf("inv....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        m = n_randint(10);
        mod = n_nextprime(n_randtest_not_zero() - 1, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);
        nmod_mat_init(C, m, m, mod);
        nmod_mat_init(I, m, m, mod);

        for (j = 0; j < m; j++)
            I->rows[j][j] = 1UL;

        /* Verify that A * A^-1 = I for random matrices */

        nmod_mat_randrank(A, m);
        /* Dense or sparse? */
        if (n_randint(2))
            nmod_mat_randops(A, 1+n_randint(1+m*m));

        result = nmod_mat_inv(B, A);
        nmod_mat_mul(C, A, B);

        if (!nmod_mat_equal(C, I) || !result)
        {
            printf("FAIL:\n");
            printf("A * A^-1 != I!\n");
            printf("A:\n");
            nmod_mat_print_pretty(A);
            printf("A^-1:\n");
            nmod_mat_print_pretty(B);
            printf("A * A^-1:\n");
            nmod_mat_print_pretty(C);
            printf("\n");
            abort();
        }

        /* Test aliasing */
        nmod_mat_set(C, A);
        nmod_mat_inv(A, A);
        nmod_mat_mul(B, A, C);

        if (!nmod_mat_equal(B, I))
        {
            printf("FAIL:\n");
            printf("aliasing failed!\n");
            nmod_mat_print_pretty(C);
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(I);
    }

    /* Test singular systems */
    for (i = 0; i < 10000; i++)
    {
        m = 1 + n_randint(10);
        r = n_randint(m);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);

        nmod_mat_randrank(A, r);

        /* Dense */
        if (n_randint(2))
            nmod_mat_randops(A, 1+n_randint(1+m*m));

        result = nmod_mat_inv(B, A);

        if (result)
        {
            printf("FAIL:\n");
            printf("singular matrix reported as invertible\n");
            abort();
        }

        /* Aliasing */
        result = nmod_mat_inv(A, A);
        if (result)
        {
            printf("FAIL:\n");
            printf("singular matrix reported as invertiblen");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    printf("PASS\n");
    return 0;
}
