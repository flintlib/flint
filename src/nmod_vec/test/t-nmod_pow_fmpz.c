/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "fmpz.h"

TEST_FUNCTION_START(nmod_vec_nmod_pow_fmpz, state)
{
    int i;

    /* check nmod_pow_fmpz matches nmod_pow_ui */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        mp_limb_t b, c1, c2;
        ulong exp1;
        fmpz_t exp2;

        nmod_init(&mod, n_randtest_not_zero(state));

        exp1 = n_randtest(state);
        fmpz_init_set_ui(exp2, exp1);

        b = n_randlimb(state) % mod.n;
        c1 = nmod_pow_ui(b, exp1, mod);
        c2 = nmod_pow_fmpz(b, exp2, mod);

        if (c1 != c2)
        {
            printf("FAIL\n");
            flint_printf("check nmod_pow_fmpz matches nmod_pow_ui\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp2);
    }

    /* check b^e1*b^e2 = b^(e1+e2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        mp_limb_t b, c1, c2, c3;
        fmpz_t exp1, exp2, exp3;

        nmod_init(&mod, n_randtest_not_zero(state));

        fmpz_init(exp1);
        fmpz_init(exp2);
        fmpz_init(exp3);

        fmpz_randtest_unsigned(exp1, state, 500);
        fmpz_randtest_unsigned(exp2, state, 500);
        fmpz_add(exp3, exp1, exp2);

        b = n_randlimb(state) % mod.n;
        c1 = nmod_pow_fmpz(b, exp1, mod);
        c2 = nmod_pow_fmpz(b, exp2, mod);
        c3 = nmod_pow_fmpz(b, exp3, mod);

        if (c3 != nmod_mul(c1, c2, mod))
        {
            printf("FAIL\n");
            flint_printf("check b^e1*b^e2 = b^(e1+e2)\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp3);
        fmpz_clear(exp2);
        fmpz_clear(exp1);
    }

    TEST_FUNCTION_END(state);
}
