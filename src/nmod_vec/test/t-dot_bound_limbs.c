/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_dot_bound_limbs, state)
{
    int i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        mp_limb_t m;
        int limbs1, limbs2;
        mpz_t t;

        len = n_randint(state, 10000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        limbs1 = _nmod_vec_dot_bound_limbs(len, mod);

        mpz_init2(t, 4*FLINT_BITS);
        flint_mpz_set_ui(t, m-1);
        mpz_mul(t, t, t);
        flint_mpz_mul_ui(t, t, len);
        limbs2 = mpz_size(t);

        if (limbs1 != limbs2)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "limbs1 = %d\n"
                    "limbs2 = %d\n"
                    "bound: %{mpz}\n",
                    m, len, limbs1, limbs2, t);

        mpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
