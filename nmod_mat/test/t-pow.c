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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main()
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("pow....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D, temp;
        mp_limb_t mod;
        slong m, j;

        do{
        m = n_randint(state, 20);
        mod=n_randint(state, 0);
        } while(m == 0 || mod == 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);
        nmod_mat_init(D, m, m, mod);
        nmod_mat_randtest(A, state);
        nmod_mat_init_set(C, A);

        nmod_mat_pow(B, A, 15);
        for(j = 1; j < 15 ; j++)
        {
            nmod_mat_mul(D, C, A);
            *temp = *C;
            *C = *D;
            *D = *temp;
        }

        if(!nmod_mat_equal(C, B))
        {
            flint_printf("FAIL: results not equal\n");fflush(stdout);
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
