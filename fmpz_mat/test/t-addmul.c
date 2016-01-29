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

    Copyright (C) 2016 Aaditya Thakkar

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);
    

    flint_printf("addmul....");
    fflush(stdout);

    for (i = 0; i < 20 *10; i++)
    {
        fmpz_mat_t A, B, C, D, T, E;

        slong m, k, n;

        m = n_randint(state, 100);
        k = n_randint(state, 100);
        n = n_randint(state, 100);

        /* Force Strassen test */
        if (i < 5)
        {
            m += 300;
            k += 300;
            n += 300;
        }

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m, n);
        fmpz_mat_init(T, m, n);
        fmpz_mat_init(E, m, n);

        fmpz_mat_randtest(A, state,n_randint(state, 200)+1);
        fmpz_mat_randtest(B, state,n_randint(state, 200)+1);
        fmpz_mat_randtest(C, state,n_randint(state, 200)+1);

        fmpz_mat_addmul(D, C, A, B);

        fmpz_mat_mul(T, A, B);
        fmpz_mat_add(E, C, T);

        if (!fmpz_mat_equal(D, E))
        {
            flint_printf("FAIL: results not equal\n");
            fmpz_mat_print_pretty(A);
            fmpz_mat_print_pretty(B);
            fmpz_mat_print_pretty(C);
            fmpz_mat_print_pretty(D);
            fmpz_mat_print_pretty(E);
            abort();
        }

        /* Check aliasing */
        fmpz_mat_addmul(C, C, A, B);

        if (!fmpz_mat_equal(C, E))
        {
            flint_printf("FAIL: results not equal (aliasing)\n");
            fmpz_mat_print_pretty(A);
            fmpz_mat_print_pretty(B);
            fmpz_mat_print_pretty(C);
            fmpz_mat_print_pretty(D);
            fmpz_mat_print_pretty(E);
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
        fmpz_mat_clear(T);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
