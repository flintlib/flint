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

    printf("poly_nmod_mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        ring_t ZZ, ZZp, ZZpx[3];
        ring_struct * ring;
        elem_ptr p;
        long size[3];
        elem_poly_t A, B, C, D, E;

        ring_init_limb(ZZ);
        ELEM_TMP_INIT(p, ZZ);
        elem_set_ui(p, n_randtest_prime(state, 0), ZZ);
        ring_init_mod(ZZp, ZZ, p);
        ring_init_poly(ZZpx[0], ZZp);
        ring_init_poly(ZZpx[1], ZZpx[0]);
        ring_init_poly(ZZpx[2], ZZpx[1]);

        ring = ZZpx[n_randint(state, 3)];

        size[0] = 1 + n_randint(state, 15);
        size[1] = 1 + n_randint(state, 15);
        size[2] = 1 + n_randint(state, 15);

        elem_init(A, ring);
        elem_init(B, ring);
        elem_init(C, ring);
        elem_init(D, ring);
        elem_init(E, ring);

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);

        elem_poly_nmod_mul(C, A, B, ring);
        elem_poly_nmod_mul(A, A, B, ring);
        if (!elem_equal(C, A, ring))
        {
            printf("FAIL: aliasing C, A\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_poly_nmod_mul(C, A, B, ring);
        elem_poly_nmod_mul(B, A, B, ring);
        if (!elem_equal(C, B, ring))
        {
            printf("FAIL: aliasing C, B\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_set(B, A, ring);
        elem_poly_nmod_mul(C, A, B, ring);
        elem_poly_nmod_mul(D, A, A, ring);
        if (!elem_equal(C, D, ring))
        {
            printf("FAIL: aliasing A, B\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_poly_nmod_mul(B, A, A, ring);
        elem_poly_nmod_mul(A, A, A, ring);
        if (!elem_equal(B, A, ring))
        {
            printf("FAIL: aliasing C, A, B\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_randtest(C, state, size, ring);
        elem_poly_nmod_mul(D, A, B, ring);
        elem_poly_nmod_mul(D, D, C, ring);
        elem_poly_nmod_mul(E, B, C, ring);
        elem_poly_nmod_mul(E, A, E, ring);
        if (!elem_equal(D, E, ring))
        {
            printf("FAIL: (A*B)*C = A*(B*C)\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            elem_print(E, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_randtest(C, state, size, ring);
        elem_poly_nmod_mul(D, A, B, ring);
        elem_poly_nmod_mul(E, A, C, ring);
        elem_add(E, D, E, ring);
        elem_add(D, B, C, ring);
        elem_poly_nmod_mul(D, A, D, ring);

        if (!elem_equal(D, E, ring))
        {
            printf("FAIL: A*(B+C) = A*B+A*C\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            elem_print(E, ring); printf("\n\n");
            abort();
        }

        elem_clear(A, ring);
        elem_clear(B, ring);
        elem_clear(C, ring);
        elem_clear(D, ring);
        elem_clear(E, ring);

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

