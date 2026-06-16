/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix_padic.h"
#include "gr.h"
#include "fmpz.h"

TEST_FUNCTION_START(radix_padic, state)
{
    slong iter;

    /* Generic ring tests over a range of primes / precision configurations. */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        ulong p = n_randtest_prime(state, 0);
        slong prec_abs, prec_rel;
        int flags = 0;

        if (n_randint(state, 2))
            flags |= RADIX_PADIC_SIGNED;

        switch (n_randint(state, 4))
        {
            case 0:   /* classic padic: absolute precision only */
                prec_abs = 1 + n_randint(state, 20);
                prec_rel = RADIX_PADIC_PREC_INF;
                break;
            case 1:   /* relative precision only */
                prec_abs = RADIX_PADIC_PREC_INF;
                prec_rel = 1 + n_randint(state, 20);
                break;
            case 2:   /* both */
                prec_abs = 1 + n_randint(state, 20);
                prec_rel = 1 + n_randint(state, 20);
                break;
            default:  /* exact ring */
                prec_abs = RADIX_PADIC_PREC_INF;
                prec_rel = RADIX_PADIC_PREC_INF;
                break;
        }

        if (n_randint(state, 2))
            flags |= RADIX_PADIC_TEST_LIMITS;

        gr_ctx_init_radix_padic(ctx, p, prec_rel, prec_abs, flags);

        /* allow get_fmpz roundtrip to return unable for huge test values */
        if (flags & RADIX_PADIC_TEST_LIMITS)
            ctx->size_limit = (1 << 16);

        gr_test_ring(ctx, 30, 0);

        /* Check optimized add vs reference implementation. */
        for (slong j = 0; j < 100; j++)
        {
            radix_padic_t x, y, r1, r2;
            int s1, s2;
            int aliasing;

            radix_padic_init(x, ctx);
            radix_padic_init(y, ctx);
            radix_padic_init(r1, ctx);
            radix_padic_init(r2, ctx);

            GR_IGNORE(gr_randtest(x, state, ctx));
            GR_IGNORE(gr_randtest(y, state, ctx));

            GR_IGNORE(gr_randtest(r1, state, ctx));
            GR_IGNORE(gr_randtest(r2, state, ctx));

            aliasing = n_randint(state, 5);

            switch (aliasing)
            {
                case 0:
                    s1 = radix_padic_add(r1, x, y, ctx);
                    break;
                case 1:
                    radix_padic_set(r1, x, ctx);
                    s1 = radix_padic_add(r1, r1, y, ctx);
                    break;
                case 2:
                    radix_padic_set(r1, y, ctx);
                    s1 = radix_padic_add(r1, x, r1, ctx);
                    break;
                case 3:
                    radix_padic_set(x, y, ctx);
                    s1 = radix_padic_add(r1, x, x, ctx);
                    break;
                default:
                    radix_padic_set(x, y, ctx);
                    radix_padic_set(r1, x, ctx);
                    s1 = radix_padic_add(r1, r1, r1, ctx);
                    break;
            }

            s2 = _radix_padic_add_sub_reference(r2, x, y, 0, ctx);

            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS &&
                (r1->v != r2->v || r1->N != r2->N ||
                !radix_integer_equal(&r1->u, &r2->u, RADIX_PADIC_CTX_RADIX(ctx))))
            {
                flint_printf("FAIL: p-adic add vs reference\n");
                flint_printf("s1 = %d\n", s1);
                flint_printf("s2 = %d\n", s2);
                flint_printf("x = "); gr_println(x, ctx);
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("r1 = "); gr_println(r1, ctx);
                flint_printf("r2 = "); gr_println(r2, ctx);
                flint_abort();
            }

            radix_padic_clear(x, ctx);
            radix_padic_clear(y, ctx);
            radix_padic_clear(r1, ctx);
            radix_padic_clear(r2, ctx);
        }

        /* Check optimized mul vs reference implementation. */
        for (slong j = 0; j < 100; j++)
        {
            radix_padic_t x, y, r1, r2;
            int s1, s2;
            int aliasing;

            radix_padic_init(x, ctx);
            radix_padic_init(y, ctx);
            radix_padic_init(r1, ctx);
            radix_padic_init(r2, ctx);

            GR_IGNORE(gr_randtest(x, state, ctx));
            GR_IGNORE(gr_randtest(y, state, ctx));

            GR_IGNORE(gr_randtest(r1, state, ctx));
            GR_IGNORE(gr_randtest(r2, state, ctx));

            aliasing = n_randint(state, 5);

            switch (aliasing)
            {
                case 0:
                    s1 = radix_padic_mul(r1, x, y, ctx);
                    break;
                case 1:
                    radix_padic_set(r1, x, ctx);
                    s1 = radix_padic_mul(r1, r1, y, ctx);
                    break;
                case 2:
                    radix_padic_set(r1, y, ctx);
                    s1 = radix_padic_mul(r1, x, r1, ctx);
                    break;
                case 3:
                    radix_padic_set(x, y, ctx);
                    s1 = radix_padic_mul(r1, x, x, ctx);
                    break;
                default:
                    radix_padic_set(x, y, ctx);
                    radix_padic_set(r1, x, ctx);
                    s1 = radix_padic_mul(r1, r1, r1, ctx);
                    break;
            }

            s2 = _radix_padic_mul_reference(r2, x, y, ctx);

            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS &&
                (r1->v != r2->v || r1->N != r2->N ||
                !radix_integer_equal(&r1->u, &r2->u, RADIX_PADIC_CTX_RADIX(ctx))))
            {
                flint_printf("FAIL: p-adic mul vs reference\n");
                flint_printf("s1 = %d\n", s1);
                flint_printf("s2 = %d\n", s2);
                flint_printf("x = "); gr_println(x, ctx);
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("r1 = "); gr_println(r1, ctx);
                flint_printf("r2 = "); gr_println(r2, ctx);
                flint_abort();
            }

            radix_padic_clear(x, ctx);
            radix_padic_clear(y, ctx);
            radix_padic_clear(r1, ctx);
            radix_padic_clear(r2, ctx);
        }

        /* Check consistency of square root */
        for (slong j = 0; j < 10; j++)
        {
            radix_padic_t x, y, r1, r2;
            int s1, s2;

            radix_padic_init(x, ctx);
            radix_padic_init(y, ctx);
            radix_padic_init(r1, ctx);
            radix_padic_init(r2, ctx);

            GR_IGNORE(gr_randtest(x, state, ctx));
            GR_IGNORE(gr_randtest(r1, state, ctx));
            GR_IGNORE(gr_randtest(r2, state, ctx));

            y->N = n_randint(state, 10);
            GR_MUST_SUCCEED(radix_padic_add(y, x, y, ctx));

            s1 = radix_padic_sqrt(r1, x, ctx);
            s2 = radix_padic_sqrt(r2, y, ctx);

            if ((s1 == GR_SUCCESS && s2 == GR_SUCCESS && radix_padic_equal(r1, r2, ctx) == T_FALSE)
                || (s1 == GR_SUCCESS && s2 == GR_DOMAIN) || (s1 == GR_DOMAIN && s2 == GR_SUCCESS))
            {
                flint_printf("FAIL: p-adic sqrt\n");
                flint_printf("s1 = %d\n", s1);
                flint_printf("s2 = %d\n", s2);
                flint_printf("x = "); gr_println(x, ctx);
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("r1 = "); gr_println(r1, ctx);
                flint_printf("r2 = "); gr_println(r2, ctx);
                flint_abort();
            }

            radix_padic_clear(x, ctx);
            radix_padic_clear(y, ctx);
            radix_padic_clear(r1, ctx);
            radix_padic_clear(r2, ctx);
        }


        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

