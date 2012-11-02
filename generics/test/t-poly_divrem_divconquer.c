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

    printf("poly_divrem_divconquer....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        ring_t ZZ, ZZp, ZZpx[3];
        ring_struct * ring;
        elem_ptr p;
        long size[3];
        elem_poly_t A, B, C, Q, Q2, R, R2;

        ring_init_limb(ZZ);
        ELEM_TMP_INIT(p, ZZ);
        elem_set_ui(p, n_randtest_prime(state, 0), ZZ);
        ring_init_mod(ZZp, ZZ, p);
        ring_init_poly(ZZpx[0], ZZp);
        ring_init_poly(ZZpx[1], ZZpx[0]);
        ring_init_poly(ZZpx[2], ZZpx[1]);

        ring = ZZpx[n_randint(state, 2)];

        size[0] = 1 + n_randint(state, 30);
        size[1] = 1 + n_randint(state, 30);
        size[2] = 1 + n_randint(state, 30);

        elem_init(A, ring);
        elem_init(B, ring);
        elem_init(C, ring);
        elem_init(Q, ring);
        elem_init(Q2, ring);
        elem_init(R, ring);
        elem_init(R2, ring);

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_mul(C, A, B, ring);

        elem_poly_divrem_divconquer(Q, R, C, B, ring);
        if (!elem_equal(Q, A, ring) || !elem_is_zero(R, ring))
        {
            printf("FAIL: (A * B) / B = A\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_divrem(Q, R, A, B, ring);
        elem_mul(C, Q, B, ring);
        elem_add(C, C, R, ring);
        if (!elem_equal(C, A, ring))
        {
            printf("FAIL: Q * B + R = A\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_poly_divrem_divconquer(Q, R, A, B, ring);
        elem_poly_divrem_divconquer(A, R2, A, B, ring);
        if (!elem_equal(A, Q, ring) || !elem_equal(R, R2, ring))
        {
            printf("FAIL: aliasing Q, A\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            elem_print(R2, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_poly_divrem_divconquer(Q, R, A, B, ring);
        elem_poly_divrem_divconquer(Q2, A, A, B, ring);
        if (!elem_equal(A, R, ring) || !elem_equal(Q, Q2, ring))
        {
            printf("FAIL: aliasing R, A\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(Q2, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_poly_divrem_divconquer(Q, R, A, B, ring);
        elem_poly_divrem_divconquer(Q2, B, A, B, ring);
        if (!elem_equal(B, R, ring) || !elem_equal(Q, Q2, ring))
        {
            printf("FAIL: aliasing R, B\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(Q2, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest_not_zero(B, state, size, ring);
        elem_poly_divrem_divconquer(Q, R, A, B, ring);
        elem_poly_divrem_divconquer(B, R2, A, B, ring);
        if (!elem_equal(B, Q, ring) || !elem_equal(R, R2, ring))
        {
            printf("FAIL: aliasing Q, B\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(Q, ring); printf("\n\n");
            elem_print(R, ring); printf("\n\n");
            elem_print(R2, ring); printf("\n\n");
            abort();
        }

        elem_clear(A, ring);
        elem_clear(B, ring);
        elem_clear(C, ring);
        elem_clear(Q, ring);
        elem_clear(Q2, ring);
        elem_clear(R, ring);
        elem_clear(R2, ring);

        ring_clear(ZZpx[2]);
        ring_clear(ZZpx[1]);
        ring_clear(ZZpx[0]);
        ring_clear(ZZp);
        ELEM_TMP_CLEAR(p, ZZ);
        ring_clear(ZZ);
    }

    printf("PASS\n");
    flint_randclear(state);
    return EXIT_SUCCESS;
}

