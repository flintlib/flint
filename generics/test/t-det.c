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

int bad_frac_type(ring_t r)
{
    if (r->type == TYPE_FRAC)
    {
        if (RING_NUMER(r) == RING_DENOM(r))
            return 0;

        if (RING_NUMER(r)->type == TYPE_POLY &&
            RING_PARENT(RING_NUMER(r)) == RING_DENOM(r))
            return 0;

        return 1;
    }

    if (r->type == TYPE_POLY)
        return bad_frac_type(RING_PARENT(r));

    return 0;
}

int main()
{
    flint_rand_t state;
    long iter;

    printf("det....");

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int i, depth;
        ring_t rings[5];
        ring_t mat_ring;
        long size[5] = {4, 4, 4, 4, 4};
        long n;
        elem_mat_t A, B, C;
        elem_ptr a, b, ab, c;

        /* XXX: we want to properly generate integral domains (?) */
        do {
            depth = ring_init_randtest(rings, state, 4);
        } while (bad_frac_type(rings[depth-1]) ||
            (rings[0]->type == TYPE_LIMB &&
            (depth == 1 ||
            (depth > 1 && rings[1]->type != TYPE_MOD))));

        ring_init_mat(mat_ring, rings[depth-1]);

        n = n_randint(state, 6);
        elem_mat_init(A, n, n, mat_ring);
        elem_mat_init(B, n, n, mat_ring);
        elem_mat_init(C, n, n, mat_ring);

        ELEM_TMP_INIT(a, rings[depth-1]);
        ELEM_TMP_INIT(b, rings[depth-1]);
        ELEM_TMP_INIT(ab, rings[depth-1]);
        ELEM_TMP_INIT(c, rings[depth-1]);

        elem_mat_randtest(A, state, size, mat_ring);
        elem_mat_randtest(B, state, size, mat_ring);

        elem_mat_mul(C, A, B, mat_ring);

        elem_mat_det(a, A, mat_ring);
        elem_mat_det(b, B, mat_ring);
        elem_mat_det(c, C, mat_ring);

        elem_mul(ab, a, b, rings[depth-1]);

        if (!elem_equal(ab, c, rings[depth-1]))
        {
            printf("FAIL: det(A)*det(B) != det(A * B)\n");
            ring_print(rings[depth-1]); printf("\n\n");
            ring_print(mat_ring); printf("\n\n");
            elem_mat_print(A, mat_ring); printf("\n\n");
            elem_mat_print(B, mat_ring); printf("\n\n");
            elem_mat_print(C, mat_ring); printf("\n\n");
            elem_print(a, rings[depth-1]); printf("\n\n");
            elem_print(b, rings[depth-1]); printf("\n\n");
            elem_print(ab, rings[depth-1]); printf("\n\n");
            elem_print(c, rings[depth-1]); printf("\n\n");
            abort();
        }

        ELEM_TMP_CLEAR(a, rings[depth-1]);
        ELEM_TMP_CLEAR(b, rings[depth-1]);
        ELEM_TMP_CLEAR(ab, rings[depth-1]);
        ELEM_TMP_CLEAR(c, rings[depth-1]);

        elem_mat_clear(A, mat_ring);
        elem_mat_clear(B, mat_ring);
        elem_mat_clear(C, mat_ring);

        for (i = 0; i < depth; i++)
            ring_clear(rings[i]);

        ring_clear(mat_ring);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

