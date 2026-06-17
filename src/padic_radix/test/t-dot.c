/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "padic_radix.h"
#include "gr.h"
#include "gr_vec.h"
#include "fmpz.h"

TEST_FUNCTION_START(padic_radix_dot, state)
{
    slong iter;

    /* Complementary to the generic tests, check high precision */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        ulong p = n_randtest_prime(state, 0);
        slong prec_abs, prec_rel, len;
        int flags = 0;
        padic_radix_struct * v1, * v2;
        padic_radix_t r1, r2;
        int s1, s2, initial, negate;
        slong i;

        if (n_randint(state, 2))
            flags |= PADIC_RADIX_SIGNED;

        prec_abs = 1 + n_randint(state, 1000);
        prec_rel = 1 + n_randint(state, 1000);

        if (n_randint(state, 2))
            flags |= PADIC_RADIX_TEST_LIMITS;

        gr_ctx_init_padic_radix(ctx, p, prec_rel, prec_abs, flags);

        len = n_randint(state, 8);

        v1 = gr_heap_init_vec(len, ctx);
        v2 = gr_heap_init_vec(len, ctx);
        padic_radix_init(r1, ctx);
        padic_radix_init(r2, ctx);

        for (i = 0; i < len; i++)
        {
            radix_integer_rand_limbs(PADIC_RADIX_UNIT(v1 + i), state, 300, PADIC_RADIX_CTX_RADIX(ctx));
            radix_integer_rand_limbs(PADIC_RADIX_UNIT(v2 + i), state, 300, PADIC_RADIX_CTX_RADIX(ctx));

            PADIC_RADIX_VAL(v1 + i) = n_randint(state, 30);
            PADIC_RADIX_VAL(v2 + i) = n_randint(state, 30);

            PADIC_RADIX_N(v1 + i) = n_randint(state, 2) ? PADIC_RADIX_EXACT : n_randint(state, 1000);
            PADIC_RADIX_N(v2 + i) = n_randint(state, 2) ? PADIC_RADIX_EXACT : n_randint(state, 1000);

            if (n_randint(state, 2))
                PADIC_RADIX_UNIT(v1 + i)->size *= -1;
            if (n_randint(state, 2))
                PADIC_RADIX_UNIT(v2 + i)->size *= -1;

            _padic_radix_finalize(v1 + i, ctx);
            _padic_radix_finalize(v2 + i, ctx);
        }

        initial = n_randint(state, 2);
        negate = n_randint(state, 2);

        padic_radix_randtest(r1, state, ctx);
        padic_radix_randtest(r2, state, ctx);

        s1 = padic_radix_dot_strided_naive(r1, initial ? r2 : NULL, negate, v1, 1, v2, 1, len, ctx);
        s2 = padic_radix_dot_strided_delayed(r2, initial ? r2 : NULL, negate, v1, 1, v2, 1, len, ctx);

        if (s1 == GR_SUCCESS && s2 == GR_SUCCESS && padic_radix_equal(r1, r2, ctx) == T_FALSE)
        {
            flint_printf("FAIL: naive vs delayed\n");
            flint_printf("r1 = "); gr_println(r1, ctx);
            flint_printf("r2 = "); gr_println(r2, ctx);
            flint_abort();
        }

        padic_radix_clear(r1, ctx);
        padic_radix_clear(r2, ctx);

        gr_heap_clear_vec(v1, len, ctx);
        gr_heap_clear_vec(v2, len, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

