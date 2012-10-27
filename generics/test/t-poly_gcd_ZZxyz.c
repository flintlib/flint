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

    printf("poly_gcd_ZZxyz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ring_t ZZ, ZZx, ZZxy, ZZxyz;
        elem_poly_t A, B, C, AC, BC, G, Q, R;
        ring_struct * ring;
        long size[5] = {7, 7, 7, 7, 7};

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_poly(ZZxy, ZZx);
        ring_init_poly(ZZxyz, ZZxy);

        switch (n_randint(state, 3))
        {
            case 0:
                ring = ZZx;
                break;
            case 1:
                ring = ZZxy;
                break;
            default:
                ring = ZZxyz;
                size[0] = 4;
                size[1] = 5;
        }

        elem_init(A, ring);
        elem_init(B, ring);
        elem_init(C, ring);
        elem_init(AC, ring);
        elem_init(BC, ring);
        elem_init(G, ring);
        elem_init(Q, ring);
        elem_init(R, ring);

        elem_randtest_not_zero(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_randtest_not_zero(C, state, size, ring);

        elem_mul(AC, A, C, ring);
        elem_mul(BC, B, C, ring);

        elem_gcd(G, AC, BC, ring);
        elem_divrem(Q, R, G, C, ring);

        if (!elem_is_zero(R, ring))
        {
            printf("FAIL:\n");
            elem_print(A, ring); printf("\n");
            elem_print(B, ring); printf("\n");
            elem_print(AC, ring); printf("\n");
            elem_print(BC, ring); printf("\n");
            elem_print(C, ring); printf("\n");
            elem_print(G, ring); printf("\n");
            elem_print(R, ring); printf("\n");
            abort();
        }

        elem_clear(A, ring);
        elem_clear(B, ring);
        elem_clear(C, ring);
        elem_clear(AC, ring);
        elem_clear(BC, ring);
        elem_clear(G, ring);
        elem_clear(Q, ring);
        elem_clear(R, ring);

        ring_clear(ZZxyz);
        ring_clear(ZZxy);
        ring_clear(ZZx);
        ring_clear(ZZ);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
