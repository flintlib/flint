/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "n_mod.h"
#include "n_mod_vec.h"

TEST_FUNCTION_START(n_mod_vec_aors, state)
{
    slong ix;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        slong len = 1 + n_randint(state, 100);
        nn_ptr rp, ap, bp;
        n_mod_ctx_t ctx;
        int type, result;

        rp = _n_mod_vec_init(len);
        ap = _n_mod_vec_init(len);
        bp = _n_mod_vec_init(len);

        n_mod_ctx_init_rand(ctx, state);

        _n_mod_vec_rand(ap, state, len, ctx);
        _n_mod_vec_rand(bp, state, len, ctx);

        type = n_randint(state, 2);

        switch (type)
        {
            case 0: /* (a + b) - b = a, left aliasing */
                _n_mod_vec_set(rp, ap, len);
                _n_mod_vec_add(ap, ap, bp, len, ctx->nu);
                _n_mod_vec_sub(ap, ap, bp, len, ctx->nu);
                result = (_n_mod_vec_equal(rp, ap, len) != 0);
                break;

            case 1: /* a + (-b) = a - b, right aliasing */
                _n_mod_vec_set(rp, bp, len);
                _n_mod_vec_neg(bp, bp, len, ctx->nu);
                _n_mod_vec_add(bp, ap, bp, len, ctx->nu);
                _n_mod_vec_sub(rp, ap, rp, len, ctx->nu);
                result = (_n_mod_vec_equal(rp, bp, len) != 0);
                break;

            default: FLINT_UNREACHABLE;
        }

        if (!result)
            TEST_FUNCTION_FAIL(
                    "ix = %wd\n"
                    "type = %d\n"
                    "len = %wd\n",
                    ix, type, len);

        _n_mod_vec_clear(rp);
        _n_mod_vec_clear(ap);
        _n_mod_vec_clear(bp);

        n_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
