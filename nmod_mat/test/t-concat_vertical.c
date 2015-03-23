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

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int main(void)
{
    nmod_mat_t A, B, C;
    nmod_mat_t window1, window2;
    slong i;
    FLINT_TEST_INIT(state);


    flint_printf("concat_vertical....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong r1, r2, c1, n;

        r1 = n_randint(state, 50);
        r2 = n_randint(state, 50);
        c1 = n_randint(state, 50);
        n = n_randint(state, 50) + 1;

        nmod_mat_init(A, r1, c1, n);
        nmod_mat_init(B, r2, c1, n);
        nmod_mat_init(C, (r1+r2), c1, n);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_randtest(C, state);

        nmod_mat_concat_vertical(C, A, B);
        
        nmod_mat_window_init(window1, C, 0, 0, r1, c1);
        nmod_mat_window_init(window2, C, r1, 0, (r1+r2), c1);

        if (!(nmod_mat_equal(window1, A) && nmod_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            nmod_mat_print_pretty(A);
            flint_printf("B = \n");
            nmod_mat_print_pretty(B);
            flint_printf("A concat_vertical B = \n");
            nmod_mat_print_pretty(C);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        
        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);

        nmod_mat_window_clear(window1);
        nmod_mat_window_clear(window2);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
