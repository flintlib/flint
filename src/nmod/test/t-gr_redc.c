/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_nmod_redc, state)
{
    gr_ctx_t ZZn;
    /* int flags = GR_TEST_ALWAYS_ABLE; */
    int flags = 0;
    slong iter;
    ulong n;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        if (n_randint(state, 2))
            n = (n_randbits(state, FLINT_BITS - 2) | 1);
        else
            n = n_randtest(state);

        if (gr_ctx_init_nmod_redc_fast(ZZn, n) == GR_SUCCESS)
        {
            if (n_randint(state, 2))
                GR_MUST_SUCCEED(gr_ctx_set_is_field(ZZn, n_is_prime(n) ? T_TRUE : T_FALSE));

            _gr_poly_test_mullow((gr_method_poly_binary_trunc_op) _gr_poly_mullow, NULL, state, 10, 50, ZZn);
            _gr_poly_test_divrem((gr_method_poly_binary_binary_op) _gr_poly_divrem, state, 10, 30, ZZn);

            gr_test_ring(ZZn, 100, flags);
            gr_ctx_clear(ZZn);
        }
    }

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        if (n_randint(state, 2))
            n = n_randlimb(state) | 1;
        else
            n = n_randtest(state);

        if (gr_ctx_init_nmod_redc(ZZn, n) == GR_SUCCESS)
        {
            if (n_randint(state, 2))
                GR_MUST_SUCCEED(gr_ctx_set_is_field(ZZn, n_is_prime(n) ? T_TRUE : T_FALSE));

            _gr_poly_test_mullow((gr_method_poly_binary_trunc_op) _gr_poly_mullow, NULL, state, 10, 50, ZZn);
            _gr_poly_test_divrem((gr_method_poly_binary_binary_op) _gr_poly_divrem, state, 10, 30, ZZn);

            gr_test_ring(ZZn, 100, flags);
            gr_ctx_clear(ZZn);
        }
    }

    TEST_FUNCTION_END(state);
}
