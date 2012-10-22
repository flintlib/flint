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

int ring_has_frac(ring_t r)
{
    if (r->type == TYPE_FRAC)
        return 1;

    if (r->type == TYPE_POLY)
        return ring_has_frac(RING_PARENT(r));

    return 0;
}

int main()
{
    flint_rand_t state;
    long iter;

    printf("divrem....");

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        int i, depth;
        ring_t rings[5];
        long size[5] = {6, 6, 6, 6, 6};
        gen_t A, B, C, Q, R;

        /* XXX: we want to properly generate integral domains (?) */
        do {
            depth = ring_init_randtest(rings, state, 5);
        } while (ring_has_frac(rings[depth - 1]) ||
            (rings[0]->type == TYPE_LIMB &&
            (depth == 1 ||
            (depth > 1 && rings[1]->type != TYPE_MOD))));

        gen_init(A, rings[depth-1]);
        gen_init(B, rings[depth-1]);
        gen_init(C, rings[depth-1]);
        gen_init(Q, rings[depth-1]);
        gen_init(R, rings[depth-1]);

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_mul(C, A, B);
        gen_divrem(Q, R, C, B);
        if (!gen_equal(Q, A) || !gen_is_zero(R))
        {
            printf("FAIL: (A * B) / B = A\n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_divrem(Q, R, A, B);
        gen_mul(C, Q, B);
        gen_add(C, C, R);
        if (!gen_equal(C, A))
        {
            printf("FAIL: Q * B + R = A\n");
            gen_print(A);
            gen_print(B);
            gen_print(C);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_clear(A);
        gen_clear(B);
        gen_clear(C);
        gen_clear(Q);
        gen_clear(R);

        for (i = 0; i < depth; i++)
            ring_clear(rings[i]);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

