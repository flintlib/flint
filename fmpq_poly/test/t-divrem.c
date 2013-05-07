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

    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("divrem....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of {q,r} and {a,b} */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t A, B;
        fmpq_poly_t a, b, q, r;

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_init(r);

        fmpq_poly_randtest(A, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(B, state, n_randint(state, 50) + 1, 200);
        fmpq_poly_set(a, A);
        fmpq_poly_set(b, B);

        fmpq_poly_divrem(q, r, a, b);
        fmpq_poly_divrem(a, b, a, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(a) ? 0 : 2;
        result = (fmpq_poly_equal(q, a)) && (fmpq_poly_equal(r, b)) && !cflags;
        if (!result)
        {
            printf("FAIL (aliasing {q,r} and {a,b}):\n\n");
            printf("A = "), fmpq_poly_debug(A), printf("\n\n");
            printf("B = "), fmpq_poly_debug(B), printf("\n\n");
            printf("q = "), fmpq_poly_debug(q), printf("\n\n");
            printf("r = "), fmpq_poly_debug(r), printf("\n\n");
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b = "), fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r);
    }

    /* Check aliasing of {q,r} and {b,a} */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q, r;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_divrem(q, r, a, b);
        fmpq_poly_divrem(b, a, a, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(q, b)) && (fmpq_poly_equal(r, a)) && !cflags;
        if (!result)
        {
            printf("FAIL (aliasing of {q,r} and {b,a}):\n\n");
            printf("q = "), fmpq_poly_debug(q), printf("\n\n");
            printf("r = "), fmpq_poly_debug(r), printf("\n\n");
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b = "), fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r);
    }

    /* check a = q b + r */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q, r, rhs;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_init(r);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_divrem(q, r, a, b);
        fmpq_poly_mul(rhs, q, b);
        fmpq_poly_add(rhs, rhs, r);

        cflags |= fmpq_poly_is_canonical(q)   ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = fmpq_poly_equal(a, rhs) && !cflags;
        if (!result)
        {
            printf("FAIL (a == q b + r):\n\n");
            printf("a       = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b       = "), fmpq_poly_debug(b), printf("\n\n");
            printf("q       = "), fmpq_poly_debug(q), printf("\n\n");
            printf("r       = "), fmpq_poly_debug(r), printf("\n\n");
            printf("q b + r = "), fmpq_poly_debug(rhs), printf("\n\n");
            printf("cflags  = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r);
        fmpq_poly_clear(rhs);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
