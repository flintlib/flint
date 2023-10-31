/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_scalar_addmul_nmod, state)
{
    int i, result;

    /* Check (a + b*c) == a + (b*c) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len = n_randint(state, 100) + 1;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t c = n_randint(state, n);
        nmod_t mod;

        mp_ptr vec = _nmod_vec_init(len);
        mp_ptr vec2 = _nmod_vec_init(len);
        mp_ptr vec3 = _nmod_vec_init(len);

        nmod_init(&mod, n);

        _nmod_vec_randtest(vec, state, len, mod);
        _nmod_vec_randtest(vec2, state, len, mod);
        flint_mpn_copyi(vec3, vec2, len);

        _nmod_vec_scalar_mul_nmod(vec3, vec, len, c, mod);
        _nmod_vec_add(vec3, vec3, vec2, len, mod);

        _nmod_vec_scalar_addmul_nmod(vec2, vec, len, c, mod);

        result = _nmod_vec_equal(vec2, vec3, len);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("len = %wd, n = %wd\n", len, n);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(vec);
        _nmod_vec_clear(vec2);
        _nmod_vec_clear(vec3);
    }

    TEST_FUNCTION_END(state);
}
