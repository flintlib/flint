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

TEST_FUNCTION_START(arith_bell_number_nmod_vec, state)
{
    mp_ptr b1, b2, b3;
    slong n;

    const slong maxn = 3000;


    b1 = _nmod_vec_init(maxn);
    b2 = _nmod_vec_init(maxn);
    b3 = _nmod_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        nmod_t mod;
        mp_limb_t p;

        p = n_randtest_not_zero(state);
        nmod_init(&mod, p);

        arith_bell_number_nmod_vec_recursive(b1, n, mod);
        arith_bell_number_nmod_vec_ogf(b2, n, mod);

        if (!_nmod_vec_equal(b1, b2, n))
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        if (arith_bell_number_nmod_vec_series(b3, n, mod))
        {
            if (!_nmod_vec_equal(b1, b3, n))
            {
                flint_printf("FAIL (2):\n");
                flint_printf("p = %wu\n", p);
                flint_printf("n = %wd\n", n);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    _nmod_vec_clear(b1);
    _nmod_vec_clear(b2);
    _nmod_vec_clear(b3);

    TEST_FUNCTION_END(state);
}
