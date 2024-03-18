/*
    Copyright (C) 2009, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "long_extras.h"

TEST_FUNCTION_START(n_powmod2_preinv, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, d, r1, r2, dinv;
        mpz_t a_m, d_m, r2_m;
        slong exp;

        mpz_init(a_m);
        mpz_init(d_m);
        mpz_init(r2_m);

        d = n_randtest_not_zero(state);
        do
        {
            a = n_randtest(state);
        } while (n_gcd(d, a) != 1);
        exp = z_randtest(state);

        dinv = n_preinvert_limb(d);
        r1 = n_powmod2_preinv(a, exp, d, dinv);

        flint_mpz_set_ui(a_m, a);
        flint_mpz_set_ui(d_m, d);
        if (exp < 0)
        {
            flint_mpz_powm_ui(r2_m, a_m, -exp, d_m);
            mpz_invert(r2_m, r2_m, d_m);
        }
        else
            flint_mpz_powm_ui(r2_m, a_m, exp, d_m);
        r2 = flint_mpz_get_ui(r2_m);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "a = %wu, exp = %wd, d = %wu, dinv = %wu\n"
                    "r1 = %wu, r2 = %wu\n",
                    a, exp, d, dinv, r1, r2);

        mpz_clear(a_m);
        mpz_clear(d_m);
        mpz_clear(r2_m);
    }

    /* check 0^0 = 1 */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong d, r, dinv;

        d = n_randtest_not_zero(state);

        dinv = n_preinvert_limb(d);
        r = n_powmod2_preinv(0, 0, d, dinv);

        result = (r == 1 || (d == 1 && r == 0));
        if (!result)
            TEST_FUNCTION_FAIL("0^0 != 1 mod %wd\n", d);
    }

    /* check 0^exp = 0 mod 1 */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong r, dinv;
        slong exp;

        exp = z_randtest(state);

        dinv = n_preinvert_limb(1);
        r = n_powmod2_preinv(0, exp, 1, dinv);

        result = (r == 0);
        if (!result)
            TEST_FUNCTION_FAIL("0^%wd != 0 mod 1\n", exp);
    }

    TEST_FUNCTION_END(state);
}
