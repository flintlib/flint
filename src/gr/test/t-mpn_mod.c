/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"

TEST_FUNCTION_START(gr_mpn_mod, state)
{
    gr_ctx_t ZZn, VecZZn, MatZZn, ZZnx;
    fmpz_t n;
    slong iter;
    /* int flags = GR_TEST_ALWAYS_ABLE; */
    int flags = 0;

    fmpz_init(n);

    /* test prime close to the supported size limit */
    {
        if (FLINT_BITS == 64)
        {
            fmpz_ui_pow_ui(n, 2, 1024);
            fmpz_sub_ui(n, n, 105);
        }
        else
        {
            fmpz_ui_pow_ui(n, 2, 512);
            fmpz_sub_ui(n, n, 569);
        }

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ZZn, n));
        gr_ctx_mpn_mod_set_primality(ZZn, T_TRUE);
        gr_test_ring(ZZn, 100 * flint_test_multiplier(), flags);
        gr_ctx_clear(ZZn);
    }

    /* something close to the two-word limit */
    {
        if (FLINT_BITS == 64)
        {
            fmpz_ui_pow_ui(n, 2, 128);
            fmpz_sub_ui(n, n, 159);
        }
        else
        {
            fmpz_ui_pow_ui(n, 2, 64);
            fmpz_sub_ui(n, n, 59);
        }

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ZZn, n));
        gr_ctx_mpn_mod_set_primality(ZZn, T_TRUE);
        gr_test_ring(ZZn, 100 * flint_test_multiplier(), flags);
        gr_ctx_clear(ZZn);
    }

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        for (;;)
        {
            fmpz_randtest_not_zero(n, state, 600);
            fmpz_abs(n, n);
            if (gr_ctx_init_mpn_mod(ZZn, n) == GR_SUCCESS)
                break;
        }

        if (n_randint(state, 2))
            gr_ctx_mpn_mod_set_primality(ZZn, fmpz_is_probabprime(n) ? T_TRUE : T_FALSE);

        gr_test_ring(ZZn, 5, flags);

        /* test matrices */
        if (n_randint(state, 10) == 0)
        {
            gr_ctx_init_matrix_ring(MatZZn, ZZn, 1 + n_randint(state, 5));
            gr_test_ring(MatZZn, 3, flags);
            gr_ctx_clear(MatZZn);
        }

        /* test vectors */
        if (n_randint(state, 10) == 0)
        {
            gr_ctx_init_vector_space_gr_vec(VecZZn, ZZn, 1 + n_randint(state, 5));
            gr_test_ring(VecZZn, 3, flags);
            gr_ctx_clear(VecZZn);
        }

        /* test polynomials */
        if (n_randint(state, 10) == 0)
        {
            gr_ctx_init_gr_poly(ZZnx, ZZn);
            gr_test_ring(ZZnx, 3, flags);
            gr_ctx_clear(ZZnx);
        }


        gr_ctx_clear(ZZn);
    }


    fmpz_clear(n);

    TEST_FUNCTION_END(state);
}
