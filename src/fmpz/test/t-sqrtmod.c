/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_sqrtmod, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random integers */
    {
        int ans;
        fmpz_t a, b, c, p;
        mp_limb_t prime;

        prime = n_randint(state, UWORD(1) << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        fmpz_randm(a, state, p);

        ans = fmpz_sqrtmod(b, a, p);

        fmpz_mul(c, b, b);
        fmpz_mod(c, c, p);

        result = (ans == 0 || fmpz_equal(a, c)) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL (random):\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random squares */
    {
        int ans;
        fmpz_t a, b, c, d, p;
        mp_limb_t prime;

        prime = n_randint(state, UWORD(1) << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        do
            fmpz_randm(b, state, p);
        while (fmpz_is_zero(b));

        fmpz_mul(a, b, b);
        fmpz_mod(a, a, p);

        /* check a special case */
        if (i == 0)
        {
            fmpz_set_str(p, "15951355998396157", 10);
            fmpz_set_str(a, "7009303413761286", 10);
        }

        ans = fmpz_sqrtmod(c, a, p);

        fmpz_mul(d, c, c);
        fmpz_mod(d, d, p);

        result = (ans && fmpz_equal(a, d)) && _fmpz_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL (squares):\n");
            flint_printf("p            = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a (= b^2)    = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b            = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c (= sqrt(a) = "), fmpz_print(c), flint_printf("\n");
            flint_printf("d (= c^2)    = "), fmpz_print(d), flint_printf("\n");
            flint_printf("ans          = %d\n", ans);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
