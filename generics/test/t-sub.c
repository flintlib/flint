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

    printf("sub....");

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        int i, depth;
        ring_t R[5];
        long size[5] = {6, 6, 6, 6, 6};
        gen_t A, B, C, D, E;

        depth = ring_init_randtest(R, state, 5);

        gen_init(A, R[depth-1]);
        gen_init(B, R[depth-1]);
        gen_init(C, R[depth-1]);
        gen_init(D, R[depth-1]);
        gen_init(E, R[depth-1]);

        gen_randtest(A, state, size);
        gen_randtest(B, state, size);
        gen_sub(C, A, B);
        gen_sub(A, A, B);
        if (!gen_equal(C, A))
        {
            printf("FAIL: aliasing C, A\n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest(B, state, size);
        gen_sub(C, A, B);
        gen_sub(B, A, B);
        if (!gen_equal(C, B))
        {
            printf("FAIL: aliasing C, B\n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            abort();
        }

        gen_randtest(A, state, size);
        gen_set(B, A);
        gen_sub(C, A, B);
        gen_sub(D, A, A);
        if (!gen_equal(C, D))
        {
            printf("FAIL: aliasing A, A\n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            gen_print(D);
            abort();
        }

        gen_randtest(A, state, size);
        gen_sub(B, A, A);
        gen_sub(A, A, A);
        if (!gen_equal(B, A))
        {
            printf("FAIL: aliasing C, A, A\n");
            gen_print(A);
            gen_print(B);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest(B, state, size);
        gen_randtest(C, state, size);
        gen_sub(D, A, B);
        gen_sub(D, D, C);
        gen_add(E, B, C);
        gen_sub(E, A, E);
        if (!gen_equal(D, E))
        {
            printf("FAIL: (A-B)-C = A-(C+B)\n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            gen_print(D);
            gen_print(E);
            abort();
        }

        gen_clear(A);
        gen_clear(B);
        gen_clear(C);
        gen_clear(D);
        gen_clear(E);

        for (i = 0; i < depth; i++)
            ring_clear(R[i]);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

