/*
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

    flint_printf("scalar_div_mpq....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpz_t r, s;
        mpq_t z;

        mpq_init(z);
        fmpz_init(r);
        fmpz_init(s);
        fmpz_randtest_not_zero(r, state, 100);
        fmpz_randtest_not_zero(s, state, 100);
        fmpz_get_mpz(mpq_numref(z), r);
        fmpz_get_mpz(mpq_denref(z), s);
        mpq_canonicalize(z);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpq(b, a, z);
        fmpq_poly_scalar_div_mpq(a, a, z);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            gmp_printf("z = %Qd\n\n", z);
            flint_printf("cflags = %wu\n\n", cflags);
            abort();
        }

        mpq_clear(z);
        fmpz_clear(r);
        fmpz_clear(s);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check that (a / n1) / n2 == a / (n1 * n2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, lhs, rhs;
        fmpz_t r, s;
        mpq_t z1, z2, z;

        fmpz_init(r);
        fmpz_init(s);
        mpq_init(z1);
        mpq_init(z2);
        mpq_init(z);

        fmpz_randtest_not_zero(r, state, 100);
        fmpz_randtest_not_zero(s, state, 100);
        fmpz_get_mpz(mpq_numref(z1), r);
        fmpz_get_mpz(mpq_denref(z1), s);
        mpq_canonicalize(z1);
        fmpz_randtest_not_zero(r, state, 100);
        fmpz_randtest_not_zero(s, state, 100);
        fmpz_get_mpz(mpq_numref(z2), r);
        fmpz_get_mpz(mpq_denref(z2), s);
        mpq_canonicalize(z2);
        mpq_mul(z, z1, z2);

        fmpq_poly_init(a);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpq(lhs, a, z1);
        fmpq_poly_scalar_div_mpq(lhs, lhs, z2);
        fmpq_poly_scalar_div_mpq(rhs, a, z);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (a / n1 / n2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            gmp_printf("z1 = %Qd\n\n", z1);
            gmp_printf("z2 = %Qd\n\n", z2);
            gmp_printf("z  = %Qd\n\n", z);
            fmpq_poly_debug(lhs), flint_printf("\n\n");
            fmpq_poly_debug(rhs), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            abort();
        }

        mpq_clear(z1);
        mpq_clear(z2);
        mpq_clear(z);
        fmpz_clear(r);
        fmpz_clear(s);
        fmpq_poly_clear(a);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    /* Check that (a + b) / n == a/n + b/n */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpz_t r, s;
        mpq_t z;

        fmpz_init(r);
        fmpz_init(s);
        mpq_init(z);

        fmpz_randtest_not_zero(r, state, 100);
        fmpz_randtest_not_zero(s, state, 100);
        fmpz_get_mpz(mpq_numref(z), r);
        fmpz_get_mpz(mpq_denref(z), s);
        mpq_canonicalize(z);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_div_mpq(lhs, a, z);
        fmpq_poly_scalar_div_mpq(rhs, b, z);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_div_mpq(lhs, lhs, z);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            flint_printf("FAIL ((a + b) / n):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            gmp_printf("z = %Qd\n\n", z);
            fmpq_poly_debug(lhs), flint_printf("\n\n");
            fmpq_poly_debug(rhs), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            abort();
        }

        mpq_clear(z);
        fmpz_clear(r);
        fmpz_clear(s);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
