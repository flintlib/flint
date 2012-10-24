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

void _set_to_lead_coeff(elem_poly_struct * A, const elem_poly_struct * B, const ring_t ring)
{
    _elem_poly_fit_length(A, 1, ring);
    elem_set(A->coeffs, INDEX(B->coeffs, B->length - 1, RING_PARENT(ring)->size), RING_PARENT(ring));
    _elem_poly_set_length(A, 1, ring);
}

void
_gen_poly_set_to_lead_coeff(gen_t A, const gen_t B)
{
    _set_to_lead_coeff(A->elem, B->elem, A->ring);
}

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

    printf("pseudo_divrem....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        int i, depth;
        ring_t rings[5];
        long size[5] = {6, 6, 6, 6, 6};
        ulong d, d2;
        gen_t A, B, C, D, Q, Q2, R, R2;

        /* XXX: we want to properly generate integral domains (?) */
        do {
            depth = ring_init_randtest(rings, state, 5);
        } while (rings[depth-1]->type != TYPE_POLY ||
            (bad_frac_type(rings[depth-1]) ||
            (rings[0]->type == TYPE_LIMB &&
            (depth == 1 ||
            (depth > 1 && rings[1]->type != TYPE_MOD)))));

        gen_init(A, rings[depth-1]);
        gen_init(B, rings[depth-1]);
        gen_init(C, rings[depth-1]);
        gen_init(D, rings[depth-1]);
        gen_init(Q, rings[depth-1]);
        gen_init(Q2, rings[depth-1]);
        gen_init(R, rings[depth-1]);
        gen_init(R2, rings[depth-1]);

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_mul(C, A, B);
        gen_pseudo_divrem(Q, R, &d, C, B);
        if (!gen_equal(Q, A) || !gen_is_zero(R) || d != 0)
        {
            printf("FAIL: (A * B) / B = A, d = 0\n");
            printf("%lu\n", d);
            gen_print(A);
            gen_print(B);
            gen_print(C);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_pseudo_divrem(Q, R, &d, A, B);

        gen_mul(C, Q, B);
        gen_add(C, C, R);
        _gen_poly_set_to_lead_coeff(D, B);
        gen_pow_ui(D, D, d);
        gen_mul(D, D, A);

        if (!gen_equal(C, D))
        {
            printf("FAIL: QB + R = lead(B)^d * A\n");
            printf("%lu\n", d);
            gen_print(A);
            gen_print(B);
            gen_print(Q);
            gen_print(R);
            gen_print(C);
            gen_print(D);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_pseudo_divrem(Q, R, &d, A, B);
        gen_pseudo_divrem(A, R2, &d2, A, B);
        if (!gen_equal(A, Q) || !gen_equal(R, R2) || d != d2)
        {
            printf("FAIL: aliasing Q, A\n");
            printf("%lu\n", d2);
            gen_print(A);
            gen_print(B);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_pseudo_divrem(Q, R, &d, A, B);
        gen_pseudo_divrem(Q2, A, &d2, A, B);
        if (!gen_equal(A, R) || !gen_equal(Q, Q2) || d != d2)
        {
            printf("FAIL: aliasing R, A\n");
            printf("%lu\n", d);
            gen_print(A);
            gen_print(B);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_pseudo_divrem(Q, R, &d, A, B);
        gen_pseudo_divrem(Q2, B, &d2, A, B);
        if (!gen_equal(B, R) || !gen_equal(Q, Q2) || d != d2)
        {
            printf("FAIL: aliasing R, B\n");
            printf("%lu\n", d);
            gen_print(A);
            gen_print(B);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_randtest(A, state, size);
        gen_randtest_not_zero(B, state, size);
        gen_pseudo_divrem(Q, R, &d, A, B);
        gen_pseudo_divrem(B, R2, &d2, A, B);
        if (!gen_equal(B, Q) || !gen_equal(R, R2) || d != d2)
        {
            printf("FAIL: aliasing Q, B\n");
            printf("%lu\n", d);
            gen_print(A);
            gen_print(B);
            gen_print(Q);
            gen_print(R);
            abort();
        }

        gen_clear(A);
        gen_clear(B);
        gen_clear(C);
        gen_clear(D);
        gen_clear(Q);
        gen_clear(Q2);
        gen_clear(R);
        gen_clear(R2);

        for (i = 0; i < depth; i++)
            ring_clear(rings[i]);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

