/*
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"

/* the product of the primes p <= n + 1 with (p - 1) | n, n even >= 2 */
static void
denom_naive(fmpz_t den, ulong n)
{
    ulong p;

    fmpz_one(den);
    for (p = 2; p <= n + 1; p = n_nextprime(p, 1))
        if (n % (p - 1) == 0)
            fmpz_mul_ui(den, den, p);
}

TEST_FUNCTION_START(arith_bernoulli_number_denom, state)
{
    fmpz_t s, t;
    slong n, iter;


    fmpz_init(s);
    fmpz_init(t);

    for (n = 0; n < 1000; n++)
    {
        arith_bernoulli_number_denom(t, n);
        fmpz_addmul_ui(s, t, n_nth_prime(n+1));
    }

    fmpz_set_str(t, "34549631155954474103407159", 10);

    if (!fmpz_equal(s, t))
    {
        flint_printf("FAIL: Hash disagrees with known value\n");
        fflush(stdout);
        flint_abort();
    }

    /* even n beyond the table against the definition; highly composite
       n among them */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        ulong m;

        if (iter % 4 == 0)
            m = 2 * n_randint(state, 50000);
        else if (iter % 4 == 1)
            m = 720720 * (1 + n_randint(state, 3));
        else
            m = 2 * (90 + n_randint(state, 1000));

        arith_bernoulli_number_denom(s, m);
        denom_naive(t, m);

        if (m >= 2 && !fmpz_equal(s, t))
        {
            flint_printf("FAIL: n = %wu\n", m);
            fmpz_print(s); flint_printf("\n");
            fmpz_print(t); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(s);
    fmpz_clear(t);

    TEST_FUNCTION_END(state);
}
