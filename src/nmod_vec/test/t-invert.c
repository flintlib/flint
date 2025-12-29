/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "test_helpers.h"

void _nmod_vec_randtest_not_zero(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    slong i, sparseness;

    if (n_randint(state, 2))
    {
        for (i = 0; i < len; i++)
            vec[i] = 1 + (n_randtest(state) % (mod.n - 1));
    }
    else
    {
        sparseness = 1 + n_randint(state, FLINT_MAX(2, len));

        for (i = 0; i < len; i++)
        {
            if (n_randint(state, sparseness))
                vec[i] = 1;
            else
                vec[i] = 1 + (n_randtest(state) % (mod.n - 1));
        }
    }
}

int check_nmod_vec_invert(slong len, ulong p, flint_rand_t state)
{
    nn_ptr vec, res1, res2, res3;
    nmod_t mod;

    ulong res = 0;

    nmod_init(&mod, p);

    vec = _nmod_vec_init(len);
    _nmod_vec_randtest_not_zero(vec, state, len, mod);

    res1 = _nmod_vec_init(len);
    _nmod_vec_invert_naive(res1, vec, len, mod);

    res2 = _nmod_vec_init(len);
    _nmod_vec_invert_generic(res2, vec, len, mod);

    if (!_nmod_vec_equal(res1, res2, len))
        res = 1;

    if (NMOD_CAN_USE_SHOUP(mod))
    {
        res3 = _nmod_vec_init(len);
        _nmod_vec_invert_shoup(res3, vec, len, mod);

        if (!_nmod_vec_equal(res3, res2, len))
            res = 2;
    }

    _nmod_vec_clear(vec);
    _nmod_vec_clear(res1);
    _nmod_vec_clear(res2);
    if (NMOD_CAN_USE_SHOUP(mod))
        _nmod_vec_clear(res3);

    return res;
}

TEST_FUNCTION_START(nmod_vec_invert, state)
{
    int i, result;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        ulong p = n_randtest_prime(state, 1);
        ulong len = 1 + n_randint(state, 1000);

        result = check_nmod_vec_invert(len, p, state);

        if (result)
            TEST_FUNCTION_FAIL("mod prime = %wu, len = %wu, ret code = %wu\n",
                    p, len, result);
    }

    TEST_FUNCTION_END(state);
}
