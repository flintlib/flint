/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2024 Vincent Neiger

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

TEST_FUNCTION_START(_nmod_vec_dot_params, state)
{
    int i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        ulong m;
        int nlimbs1, nlimbs2, nlimbs3;
        dot_params_t params;
        mpz_t t;

        len = n_randint(state, 10000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        params = _nmod_vec_dot_params(len, mod);
        nlimbs1 = _nmod_vec_dot_bound_limbs_from_params(len, mod, params);
        nlimbs2 = _nmod_vec_dot_bound_limbs(len, mod);

        mpz_init2(t, 4*FLINT_BITS);
        flint_mpz_set_ui(t, m-1);
        mpz_mul(t, t, t);
        flint_mpz_mul_ui(t, t, len);
        nlimbs3 = mpz_size(t);

        if (nlimbs1 != nlimbs3)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "nlimbs1(from params) = %d\n"
                    "nlimbs3(mpz) = %d\n"
                    "bound: %{mpz}\n",
                    m, len, nlimbs1, nlimbs3, t);

        if (nlimbs2 != nlimbs3)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "nlimbs2(from len+mod) = %d\n"
                    "nlimbs3(mpz) = %d\n"
                    "bound: %{mpz}\n",
                    m, len, nlimbs2, nlimbs3, t);

        mpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
