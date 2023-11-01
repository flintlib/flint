/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_cbrt_newton_iteration, state)
{
    int i, result;

    /* random n */

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);

        n = n_randtest(state);
        val = n_cbrt_newton_iteration(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);

        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans);
            fflush(stdout);
            flint_abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }
    /* type n^3 + k */

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans, bits;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);

        bits = n_randint(state, FLINT_BITS/3 + 1);
        n = n_randtest_bits(state, bits);
        n = n*n*n;
        n += (n_randint(state, 100) - 50);
        val = n_cbrt_newton_iteration(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);

        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans);
            fflush(stdout);
            flint_abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }

    /* type n^3 + 1 */

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans, bits;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);

        bits = n_randint(state, FLINT_BITS/3 + 1);
        n = n_randtest_bits(state, bits);
        n = n*n*n;
        n += 1;
        val = n_cbrt_newton_iteration(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);

        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans);
            fflush(stdout);
            flint_abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }

    /* type n^3 - 1 */

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans, bits;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);

        bits = n_randint(state, FLINT_BITS/3 + 1);
        n = n_randtest_bits(state, bits);
        n = n*n*n;
        n -= 1;
        val = n_cbrt_newton_iteration(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);

        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans);
            fflush(stdout);
            flint_abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }

    TEST_FUNCTION_END(state);
}
