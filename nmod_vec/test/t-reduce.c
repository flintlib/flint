/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
    

    flint_printf("reduce....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong j, len = n_randint(state, 100) + 1;
        mp_ptr vec = _nmod_vec_init(len);
        mp_ptr vec2 = _nmod_vec_init(len);

        mp_limb_t n = n_randtest_not_zero(state);
        nmod_t mod;
        nmod_init(&mod, n);

        for (j = 0; j < len; j++)
        {
            vec[j] = n_randtest(state);
            vec2[j] = vec[j];
        }

        _nmod_vec_reduce(vec, vec, len, mod);
        for (j = 0; j < len; j++)
            vec2[j] = n_mod2_preinv(vec2[j], mod.n, mod.ninv);

        result = _nmod_vec_equal(vec, vec2, len);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("len = %wd, n = %wd\n", len, n);
            abort();
        }

        _nmod_vec_clear(vec);
        _nmod_vec_clear(vec2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
