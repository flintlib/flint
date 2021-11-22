/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    ulong cflags = UWORD(0);

    FLINT_TEST_INIT(state);

    flint_printf("scalar_div_mpz....");
    fflush(stdout);

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
            flint_printf("FAIL (aliasing):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fmpz_print(n);
            fflush(stdout);
            flint_abort();
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
        slong n;

        n = z_randtest_not_zero(state);
        mpz_init(n1);
        flint_mpz_set_si(n1, n);

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
            flint_printf("FAIL (comparison with _si):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("%wd", n), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
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
            flint_printf("FAIL (a / n1 / n2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpz_print(n1), flint_printf("\n\n");
            fmpz_print(n2), flint_printf("\n\n");
            fmpq_poly_debug(lhs), flint_printf("\n\n");
            fmpq_poly_debug(rhs), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
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
            flint_printf("FAIL ((a + b) / n):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpz_print(n), flint_printf("\n\n");
            fmpq_poly_debug(lhs), flint_printf("\n\n");
            fmpq_poly_debug(rhs), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        mpz_clear(m);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
