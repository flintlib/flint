/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_dot_ptr, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        ulong m, res, res2;
        nn_ptr x, y;
        nn_ptr * z;
        slong j, offset;

        len = n_randint(state, 1000) + 1;
        m = n_randtest_not_zero(state);
        offset = n_randint(state, 10);

        nmod_init(&mod, m);

        x = _nmod_vec_init(len);
        y = _nmod_vec_init(len);
        z = flint_malloc(sizeof(nn_ptr) * len);

        _nmod_vec_randtest(x, state, len, mod);
        _nmod_vec_randtest(y, state, len, mod);

        for (j = 0; j < len; j++)
            z[j] = &y[j] + offset;

        const dot_params_t params = _nmod_vec_dot_params(len, mod);

        res = _nmod_vec_dot_ptr(x, z, -offset, len, mod, params);
        res2 = _nmod_vec_dot(x, y, len, mod, params);

        if (res != res2)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "method = %d\n",
                    m, len, params.method);

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
        flint_free(z);
    }

    TEST_FUNCTION_END(state);
}
