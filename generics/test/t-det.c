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

int main()
{
    flint_rand_t state;
    long iter;

    printf("det....");

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ring_t ZZ, ZZx, MM;
        ring_t mat_ring;
        long size[5] = {4, 4, 4, 4, 4};
        long n;
        elem_mat_t A, B, C;
        elem_ptr a, b, ab, c;

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_mat(MM, ZZx);

        ring_init_mat(mat_ring, MM);

        n = n_randint(state, 8);
        elem_mat_init(A, n, n, MM);
        elem_mat_init(B, n, n, MM);
        elem_mat_init(C, n, n, MM);

        ELEM_TMP_INIT(a, ZZx);
        ELEM_TMP_INIT(b, ZZx);
        ELEM_TMP_INIT(ab, ZZx);
        ELEM_TMP_INIT(c, ZZx);

        elem_mat_randtest(A, state, size, MM);
        elem_mat_randtest(B, state, size, MM);

        elem_mat_mul(C, A, B, MM);

        elem_mat_det(a, A, MM);
        elem_mat_det(b, B, MM);
        elem_mat_det(c, C, MM);

        elem_mul(ab, a, b, ZZx);

        if (!elem_equal(ab, c, ZZx))
        {
            printf("FAIL: det(A)*det(B) != det(A * B)\n");
            elem_mat_print(A, MM); printf("\n\n");
            elem_mat_print(B, MM); printf("\n\n");
            elem_mat_print(C, MM); printf("\n\n");
            elem_print(a, ZZx); printf("\n\n");
            elem_print(b, ZZx); printf("\n\n");
            elem_print(ab, ZZx); printf("\n\n");
            elem_print(c, ZZx); printf("\n\n");
            abort();
        }

        ELEM_TMP_CLEAR(a, ZZx);
        ELEM_TMP_CLEAR(b, ZZx);
        ELEM_TMP_CLEAR(ab, ZZx);
        ELEM_TMP_CLEAR(c, ZZx);

        elem_mat_clear(A, MM);
        elem_mat_clear(B, MM);
        elem_mat_clear(C, MM);

        ring_clear(MM);
        ring_clear(ZZx);
        ring_clear(ZZ);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

