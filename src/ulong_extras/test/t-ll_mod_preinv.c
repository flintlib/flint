/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(n_ll_mod_preinv, state)
{
    int i, result;
    fmpz_t n;

    /* check (m*d + r1) % d == r1 */
    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, nh, nl, r1, r2, m;

        d = n_randtest_not_zero(state);
        m = n_randtest(state);
        r1 = n_randtest(state) % d;

        /* <nh, nl> = m*d + r1 < (m + 1)*d <= 2^FLINT_BITS * d */

        umul_ppmm(nh, nl, m, d);
        add_ssaaaa(nh, nl, nh, nl, UWORD(0), r1);

        dinv = n_preinvert_limb(d);

        r2 = n_ll_mod_preinv(nh, nl, d, dinv);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "<nh, nl> = (m*d + r1) but <nh, nl> %% d != r1\n"
                    "nh = %wu, nl = %wu, d = %wu, dinv = %wu\n"
                    "r1 = %wu, r2 = %wu\n",
                    nh, nl, d, dinv, r1, r2);
    }

    fmpz_init(n);

    /* compare with fmpz_fdiv_ui */
    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, nh, nl, r1, r2;

        d = n_randtest_not_zero(state);
        nh = n_randtest(state);
        nl = n_randtest(state);

        dinv = n_preinvert_limb(d);

        /* n = <nh, nl> */
        fmpz_set_ui(n, nh);
        fmpz_mul_2exp(n, n, FLINT_BITS);
        fmpz_add_ui(n, n, nl);

        r1 = n_ll_mod_preinv(nh, nl, d, dinv);

        r2 = fmpz_fdiv_ui(n, d);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = <nh, nl> but n % d does not agree with fmpz_fdiv_ui\n"
                    "nh = %wu, nl = %wu, d = %wu, dinv = %wu\n"
                    "r1 = %wu, r2 = %wu\n",
                    nh, nl, d, dinv, r1, r2);
    }

    fmpz_clear(n);

    TEST_FUNCTION_END(state);
}
