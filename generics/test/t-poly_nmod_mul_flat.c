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

    printf("poly_nmod_mul_flat....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        ring_t ZZ, ZZp, ZZpx, ZZpxy;
        elem_ptr p;
        long size[2];
        elem_poly_t A, B, C, D, E;

        ring_init_limb(ZZ);
        ELEM_TMP_INIT(p, ZZ);
        elem_set_ui(p, n_randtest_prime(state, 0), ZZ);
        ring_init_mod(ZZp, ZZ, p);
        ring_init_poly(ZZpx, ZZp);
        ring_init_poly(ZZpxy, ZZpx);

        size[0] = 1 + n_randint(state, 20);
        size[1] = 1 + n_randint(state, 20);

        elem_init(A, ZZpxy);
        elem_init(B, ZZpxy);
        elem_init(C, ZZpxy);
        elem_init(D, ZZpxy);
        elem_init(E, ZZpxy);

        elem_randtest(A, state, size, ZZpxy);
        elem_randtest(B, state, size, ZZpxy);

        elem_poly_nmod_mul_flat(C, A, B, ZZpxy);
        elem_poly_nmod_mul_flat(A, A, B, ZZpxy);
        if (!elem_equal(C, A, ZZpxy))
        {
            printf("FAIL: aliasing C, A\n");
            elem_print(A, ZZpxy); printf("\n\n");
            elem_print(B, ZZpxy); printf("\n\n");
            elem_print(C, ZZpxy); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ZZpxy);
        elem_randtest(B, state, size, ZZpxy);
        elem_poly_nmod_mul_flat(C, A, B, ZZpxy);
        elem_poly_nmod_mul_flat(B, A, B, ZZpxy);
        if (!elem_equal(C, B, ZZpxy))
        {
            printf("FAIL: aliasing C, B\n");
            elem_print(A, ZZpxy); printf("\n\n");
            elem_print(B, ZZpxy); printf("\n\n");
            elem_print(C, ZZpxy); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ZZpxy);
        elem_set(B, A, ZZpxy);
        elem_poly_nmod_mul_flat(C, A, B, ZZpxy);
        elem_poly_nmod_mul_flat(D, A, A, ZZpxy);
        if (!elem_equal(C, D, ZZpxy))
        {
            printf("FAIL: aliasing A, B\n");
            elem_print(A, ZZpxy); printf("\n\n");
            elem_print(B, ZZpxy); printf("\n\n");
            elem_print(C, ZZpxy); printf("\n\n");
            elem_print(D, ZZpxy); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ZZpxy);
        elem_poly_nmod_mul_flat(B, A, A, ZZpxy);
        elem_poly_nmod_mul_flat(A, A, A, ZZpxy);
        if (!elem_equal(B, A, ZZpxy))
        {
            printf("FAIL: aliasing C, A, B\n");
            elem_print(A, ZZpxy); printf("\n\n");
            elem_print(B, ZZpxy); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ZZpxy);
        elem_randtest(B, state, size, ZZpxy);
        elem_randtest(C, state, size, ZZpxy);
        elem_poly_nmod_mul_flat(D, A, B, ZZpxy);
        elem_poly_nmod_mul_flat(D, D, C, ZZpxy);
        elem_poly_nmod_mul_flat(E, B, C, ZZpxy);
        elem_poly_nmod_mul_flat(E, A, E, ZZpxy);
        if (!elem_equal(D, E, ZZpxy))
        {
            printf("FAIL: (A*B)*C = A*(B*C)\n");
            elem_print(A, ZZpxy); printf("\n\n");
            elem_print(B, ZZpxy); printf("\n\n");
            elem_print(C, ZZpxy); printf("\n\n");
            elem_print(D, ZZpxy); printf("\n\n");
            elem_print(E, ZZpxy); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ZZpxy);
        elem_randtest(B, state, size, ZZpxy);
        elem_randtest(C, state, size, ZZpxy);
        elem_poly_nmod_mul_flat(D, A, B, ZZpxy);
        elem_poly_nmod_mul_flat(E, A, C, ZZpxy);
        elem_add(E, D, E, ZZpxy);
        elem_add(D, B, C, ZZpxy);
        elem_poly_nmod_mul_flat(D, A, D, ZZpxy);

        if (!elem_equal(D, E, ZZpxy))
        {
            printf("FAIL: A*(B+C) = A*B+A*C\n");
            elem_print(A, ZZpxy); printf("\n\n");
            elem_print(B, ZZpxy); printf("\n\n");
            elem_print(C, ZZpxy); printf("\n\n");
            elem_print(D, ZZpxy); printf("\n\n");
            elem_print(E, ZZpxy); printf("\n\n");
            abort();
        }

        elem_clear(A, ZZpxy);
        elem_clear(B, ZZpxy);
        elem_clear(C, ZZpxy);
        elem_clear(D, ZZpxy);
        elem_clear(E, ZZpxy);

        ring_clear(ZZpxy);
        ring_clear(ZZpx);
        ring_clear(ZZp);
        ELEM_TMP_CLEAR(p, ZZ);
        ring_clear(ZZ);
    }

    printf("PASS\n");
    flint_randclear(state);
    return EXIT_SUCCESS;
}

