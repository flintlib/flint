/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default.h"

TEST_FUNCTION_START(fq_default_inlines, state)
{
    slong i;

    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_t a, b, c;
        fmpz_t p, p1, o, o1;
        slong d;

        fmpz_init(p);
        fmpz_init(p1);
        fmpz_init(o);
        fmpz_init(o1);
        fmpz_randprime(p, state, n_randint(state, 2) ? 60 : 120, 1);
        d = n_randint(state, 4) + 1,
        fq_default_ctx_init(ctx, p, d, "a");
        fq_default_init(a, ctx);
        fq_default_init(b, ctx);
        fq_default_init(c, ctx);

        FLINT_TEST(fq_default_ctx_degree(ctx) == d);
        fq_default_ctx_prime(p1, ctx);
        FLINT_TEST(fmpz_equal(p1, p));
        fq_default_ctx_order(o, ctx);
        fmpz_pow_ui(o1, p, d);
        FLINT_TEST(fmpz_equal(o1, o));

        fq_default_randtest(a, state, ctx);
        fq_default_randtest(b, state, ctx);
        fq_default_randtest(c, state, ctx);

        FLINT_TEST(!!fq_default_is_zero(a, ctx) == !fq_default_is_invertible(a, ctx));

        fq_default_add(a, b, c, ctx);
        fq_default_sub(a, a, c, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_sub_one(a, c, ctx);
        fq_default_one(b, ctx);
        fq_default_sub(b, c, b, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_randtest(a, state, ctx);
        fq_default_randtest(b, state, ctx);
        fq_default_randtest_not_zero(c, state, ctx);

        fq_default_inv(a, c, ctx);
        fq_default_mul(a, a, c, ctx);
        FLINT_TEST(fq_default_is_one(a, ctx));

        fq_default_div(a, b, c, ctx);
        fq_default_mul(a, a, c, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_neg(a, c, ctx);
        fq_default_mul_si(b, c, -1, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_add(a, c, c, ctx);
        fq_default_mul_ui(b, c, 2, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_mul_fmpz(a, c, p, ctx);
        FLINT_TEST(fq_default_is_zero(a, ctx));

        fq_default_pow_ui(c, b, 2, ctx);
        fq_default_sqr(a, b, ctx);
        FLINT_TEST(fq_default_equal(a, c, ctx));

        fq_default_pow(a, b, o, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_mul(a, b, b, ctx);
        fq_default_sqr(c, b, ctx);
        FLINT_TEST(fq_default_equal(a, c, ctx));

        fq_default_randtest(a, state, ctx);
        fq_default_randtest(b, state, ctx);
        fq_default_randtest(c, state, ctx);

        fq_default_pth_root(a, b, ctx);
        fq_default_pow(a, a, p, ctx);
        FLINT_TEST(fq_default_equal(a, b, ctx));

        fq_default_sqr(a, b, ctx);
        FLINT_TEST(fq_default_is_square(a, ctx));
        FLINT_TEST(fq_default_sqrt(c, a, ctx));
        fq_default_sqr(c, c, ctx);
        FLINT_TEST(fq_default_equal(c, a, ctx));

        fq_default_clear(a, ctx);
        fq_default_clear(b, ctx);
        fq_default_clear(c, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
        fmpz_clear(p1);
        fmpz_clear(o);
        fmpz_clear(o1);
    }

    TEST_FUNCTION_END(state);
}
