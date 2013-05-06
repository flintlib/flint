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
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    ulong cflags = 0UL;

    printf("scalar_div_mpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpz_t n;
        mpz_t m;

        fmpz_init(n);
        mpz_init(m);
        fmpz_randtest_not_zero(n, state, 200);
        fmpz_get_mpz(m, n);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpz(b, a, m);
        fmpq_poly_scalar_div_mpz(a, a, m);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            printf("FAIL (aliasing):\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            fmpz_print(n);
            abort();
        }

        fmpz_clear(n);
        mpz_clear(m);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Compare with fmpq_poly_scalar_mul_si */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        mpz_t n1;
        long n;

        n = z_randtest_not_zero(state);
        mpz_init(n1);
        mpz_set_si(n1, n);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpz(b, a, n1);
        fmpq_poly_scalar_div_si(c, a, n);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(b, c) && !cflags);
        if (!result)
        {
            printf("FAIL (comparison with _si):\n");
            fmpq_poly_debug(a), printf("\n\n");
            printf("%ld", n), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            fmpq_poly_debug(c), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        mpz_clear(n1);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check that (a / n1) / n2 == a / (n1 * n2) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, lhs, rhs;
        fmpz_t n1, n2;
        mpz_t m1, m2, m;

        fmpz_init(n1);
        fmpz_init(n2);
        mpz_init(m1);
        mpz_init(m2);
        mpz_init(m);

        fmpz_randtest_not_zero(n1, state, 100);
        fmpz_randtest_not_zero(n2, state, 100);
        fmpz_get_mpz(m1, n1);
        fmpz_get_mpz(m2, n2);
        mpz_mul(m, m1, m2);

        fmpq_poly_init(a);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpz(lhs, a, m1);
        fmpq_poly_scalar_div_mpz(lhs, lhs, m2);
        fmpq_poly_scalar_div_mpz(rhs, a, m);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            printf("FAIL (a / n1 / n2):\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpz_print(n1), printf("\n\n");
            fmpz_print(n2), printf("\n\n");
            fmpq_poly_debug(lhs), printf("\n\n");
            fmpq_poly_debug(rhs), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpz_clear(n1);
        fmpz_clear(n2);
        mpz_clear(m1);
        mpz_clear(m2);
        mpz_clear(m);
        fmpq_poly_clear(a);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    /* Check that (a + b) / n == a/n + b/n */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpz_t n;
        mpz_t m;

        fmpz_init(n);
        mpz_init(m);

        fmpz_randtest_not_zero(n, state, 100);
        fmpz_get_mpz(m, n);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpz(lhs, a, m);
        fmpq_poly_scalar_div_mpz(rhs, b, m);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_div_mpz(lhs, lhs, m);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            printf("FAIL ((a + b) / n):\n");
            fmpq_poly_debug(a), printf("\n\n");
            fmpq_poly_debug(b), printf("\n\n");
            fmpz_print(n), printf("\n\n");
            fmpq_poly_debug(lhs), printf("\n\n");
            fmpq_poly_debug(rhs), printf("\n\n");
            printf("cflags = %lu\n\n", cflags);
            abort();
        }

        fmpz_clear(n);
        mpz_clear(m);
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
