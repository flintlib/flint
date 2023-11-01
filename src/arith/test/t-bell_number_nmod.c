/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "arith.h"

TEST_FUNCTION_START(arith_bell_number_nmod, state)
{
    slong i, j, iter;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mp_ptr b;
        slong n;
        nmod_t mod;
        mp_limb_t p, u;

        n = n_randint(state, 800);
        if (n_randint(state, 2))
            p = n_randtest_not_zero(state);
        else
            p = n_randtest_prime(state, 0);

        nmod_init(&mod, p);

        b = _nmod_vec_init(n + 1);
        arith_bell_number_nmod_vec(b, n + 1, mod);

        for (iter = 0; iter < 5; iter++)
        {
            j = n_randint(state, n + 1);
            u = arith_bell_number_nmod(j, mod);

            if (u != b[j])
            {
                flint_printf("FAIL: p = %wu, i = %wd\n", p, j);
                fflush(stdout);
                flint_abort();
            }
        }

        _nmod_vec_clear(b);
    }

    TEST_FUNCTION_END(state);
}
