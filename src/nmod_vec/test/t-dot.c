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

TEST_FUNCTION_START(nmod_vec_dot, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        ulong m, res;
        nn_ptr x, y;
        mpz_t s, t;
        slong j;

        len = n_randint(state, 1000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        x = _nmod_vec_init(len);
        y = _nmod_vec_init(len);

        _nmod_vec_randtest(x, state, len, mod);
        _nmod_vec_randtest(y, state, len, mod);

        const dot_params_t params = _nmod_vec_dot_params(len, mod);

        res = _nmod_vec_dot(x, y, len, mod, params);

        mpz_init(s);
        mpz_init(t);

        for (j = 0; j < len; j++)
        {
            flint_mpz_set_ui(t, x[j]);
            flint_mpz_addmul_ui(s, t, y[j]);
        }

        flint_mpz_mod_ui(s, s, m);

        if (flint_mpz_get_ui(s) != res)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "limbs1 = %d\n",
                    m, len, params);

        mpz_clear(s);
        mpz_clear(t);

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
    }

    TEST_FUNCTION_END(state);
}
