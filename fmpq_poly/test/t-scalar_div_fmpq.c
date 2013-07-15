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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    ulong cflags = 0UL;
    flint_rand_t state;

    printf("scalar_div_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpq_t z;

        fmpq_init(z);
        fmpq_randtest_not_zero(z, state, 100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_fmpq(b, a, z);
        fmpq_poly_scalar_div_fmpq(a, a, z);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            fmpq_print(z), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_clear(z);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check that (a / n1) / n2 == a / (n1 * n2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, lhs, rhs;
        fmpq_t z1, z2, z;

        fmpq_init(z1);
        fmpq_init(z2);
        fmpq_init(z);

        fmpq_randtest_not_zero(z1, state, 100);
        fmpq_randtest_not_zero(z2, state, 100);
        fmpq_mul(z, z1, z2);

        fmpq_poly_init(a);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_fmpq(lhs, a, z1);
        fmpq_poly_scalar_div_fmpq(lhs, lhs, z2);
        fmpq_poly_scalar_div_fmpq(rhs, a, z);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            printf("FAIL (a / n1 / n2):\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_print(z1), printf("\n\n");
            fmpq_print(z2), printf("\n\n");
            fmpq_print(z), printf("\n\n");
            fmpq_poly_debug(lhs), printf("\n\n");
            fmpq_poly_debug(rhs), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_clear(z1);
        fmpq_clear(z2);
        fmpq_clear(z);
        fmpq_poly_clear(a);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    /* Check that (a + b) / n == a/n + b/n */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpq_t z;

        fmpq_init(z);
        fmpq_randtest_not_zero(z, state, 100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_fmpq(lhs, a, z);
        fmpq_poly_scalar_div_fmpq(rhs, b, z);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_div_fmpq(lhs, lhs, z);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            printf("FAIL ((a + b) / n):\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            fmpq_print(z), printf("\n\n");
            fmpq_poly_debug(lhs), printf("\n\n");
            fmpq_poly_debug(rhs), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpq_clear(z);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
