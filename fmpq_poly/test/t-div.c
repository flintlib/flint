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

    printf("div....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of q and a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_div(q, a, b);
        fmpq_poly_div(a, a, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(a) ? 0 : 2;
        result = (fmpq_poly_equal(q, a) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("q = "), fmpq_poly_debug(q), printf("\n\n");
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b = "), fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
    }

    /* Check aliasing of q and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_div(q, a, b);
        fmpq_poly_div(b, a, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(q, b) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("q = "), fmpq_poly_debug(q), printf("\n\n");
            printf("a = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b = "), fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
    }

    /* Compare with divrem */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q, q2, r;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_init(q2);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_divrem(q, r, a, b);
        fmpq_poly_div(q2, a, b);

        cflags |= fmpq_poly_is_canonical(q)  ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(q2) ? 0 : 2;
        result = (fmpq_poly_equal(q, q2) && !cflags);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a  = "), fmpq_poly_debug(a), printf("\n\n");
            printf("b  = "), fmpq_poly_debug(b), printf("\n\n");
            printf("q  = "), fmpq_poly_debug(q), printf("\n\n");
            printf("r  = "), fmpq_poly_debug(r), printf("\n\n");
            printf("q2 = "), fmpq_poly_debug(q2), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
        fmpq_poly_clear(q2);
        fmpq_poly_clear(r);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
