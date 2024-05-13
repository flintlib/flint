/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_scalar_mul_nmod, state)
{
    int i, result;

    {
        nmod_t mod;
        nmod_init(&mod, 2);
        _nmod_vec_scalar_mul_nmod(NULL, NULL, 0, 1, mod);
    }

    /* Check (a + b)*c == a*c + b*c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len = n_randint(state, 100) + 1;
        ulong n = n_randtest_not_zero(state);
        ulong c = n_randint(state, n);
        nmod_t mod;

        nn_ptr vec = _nmod_vec_init(len);
        nn_ptr vec2 = _nmod_vec_init(len);
        nn_ptr vec3 = _nmod_vec_init(len);

        nmod_init(&mod, n);

        _nmod_vec_randtest(vec, state, len, mod);
        _nmod_vec_randtest(vec2, state, len, mod);

        _nmod_vec_add(vec3, vec, vec2, len, mod);
        _nmod_vec_scalar_mul_nmod(vec3, vec3, len, c, mod);

        _nmod_vec_scalar_mul_nmod(vec, vec, len, c, mod);
        _nmod_vec_scalar_mul_nmod(vec2, vec2, len, c, mod);
        _nmod_vec_add(vec, vec, vec2, len, mod);

        result = _nmod_vec_equal(vec, vec3, len);
        if (!result)
            TEST_FUNCTION_FAIL("len = %wd, n = %wd\n", len, n);

        _nmod_vec_clear(vec);
        _nmod_vec_clear(vec2);
        _nmod_vec_clear(vec3);
    }

    TEST_FUNCTION_END(state);
}
