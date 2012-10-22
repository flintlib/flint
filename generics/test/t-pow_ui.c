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

    printf("pow_ui....");

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        int i, depth;
        ring_t R[5];
        long size[5] = {5, 5, 5, 5, 5};
        gen_t A, B, C, D, E, F;
        ulong e;

        depth = ring_init_randtest(R, state, 4);

        if (depth < 3)
            e = n_randint(state, 9);
        else
            e = n_randint(state, 6);

        gen_init(A, R[depth-1]);
        gen_init(B, R[depth-1]);
        gen_init(C, R[depth-1]);
        gen_init(D, R[depth-1]);
        gen_init(E, R[depth-1]);
        gen_init(F, R[depth-1]);

        gen_randtest(A, state, size);
        gen_randtest(B, state, size);
        gen_pow_ui(B, A, e);
        gen_pow_ui(A, A, e);
        if (!gen_equal(B, A))
        {
            printf("FAIL: aliasing B, A\n");
            printf("e = %lu\n", e);
            gen_print(A);
            gen_print(B);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest(B, state, size);
        gen_mul(C, A, B);
        gen_pow_ui(D, C, e);
        gen_pow_ui(E, A, e);
        gen_pow_ui(F, B, e);
        gen_mul(F, E, F);

        if (!gen_equal(D, F))
        {
            printf("FAIL: (A * B)^e = A^e * B^e \n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            gen_print(D);
            gen_print(E);
            gen_print(F);
            abort();
        }

        gen_clear(A);
        gen_clear(B);
        gen_clear(C);
        gen_clear(D);
        gen_clear(E);
        gen_clear(F);

        for (i = 0; i < depth; i++)
            ring_clear(R[i]);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
