/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default_poly.h"

TEST_FUNCTION_START(fq_default_poly_inlines, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg;
        fq_default_ctx_t ctx;
        fq_default_t t;
        fq_default_poly_t a, b, c, d, e;

        fmpz_init(p);
        fmpz_randprime(p, state, n_randint(state, 2) ? 60 : 120, 1);
        deg = n_randint(state, 4) + 1;
        fq_default_ctx_init(ctx, p, deg, "a");

        fq_default_init(t, ctx);
        fq_default_poly_init(a, ctx);
        fq_default_poly_init(b, ctx);
        fq_default_poly_init(c, ctx);
        fq_default_poly_init(d, ctx);
        fq_default_poly_init(e, ctx);

        fq_default_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_default_poly_randtest(b, state, n_randint(state, 100), ctx);
        fq_default_poly_randtest(c, state, n_randint(state, 100), ctx);

        fq_default_poly_get_coeff(t, a, 0, ctx);
        FLINT_TEST(!fq_default_poly_is_unit(a, ctx) ==
                                       !fq_default_is_invertible(t, ctx) ||
                                        fq_default_poly_length(a, ctx) != 1);

        FLINT_TEST(fq_default_poly_hamming_weight(a, ctx) <=
                                               fq_default_poly_length(a, ctx));

        fq_default_poly_add(a, b, c, ctx);
        fq_default_poly_sub(a, a, c, ctx);
        FLINT_TEST(fq_default_poly_equal(a, b, ctx));

        fq_default_poly_neg(b, a, ctx);
        fq_default_poly_add(b, a, b, ctx);
        fq_default_poly_zero(a, ctx);
        FLINT_TEST(fq_default_poly_equal(a, b, ctx));

        fq_default_poly_randtest(a, state, n_randint(state, 100), ctx);
        fq_default_poly_randtest(b, state, n_randint(state, 100), ctx);
        fq_default_poly_randtest(c, state, n_randint(state, 100), ctx);
        fq_default_poly_randtest(d, state, n_randint(state, 100), ctx);
        fq_default_randtest(t, state, ctx);

        fq_default_poly_add_si(a, b, 1, ctx);
        fq_default_poly_add_si(a, a, -1, ctx);
        FLINT_TEST(fq_default_poly_equal(a, b, ctx));

        fq_default_poly_set_fq_default(d, t, ctx);
        fq_default_poly_scalar_mul_fq_default(a, c, t, ctx);
        fq_default_poly_mul(b, c, d, ctx);
        FLINT_TEST(fq_default_poly_equal(a, b, ctx));
        if (!fq_default_is_zero(t, ctx))
        {
            fq_default_poly_scalar_div_fq_default(a, a, t, ctx);
            FLINT_TEST(fq_default_poly_equal(a, c, ctx));
        }

        fq_default_poly_set(a, b, ctx);
        fq_default_poly_scalar_addmul_fq_default(a, c, t, ctx);
        fq_default_poly_scalar_mul_fq_default(d, c, t, ctx);
        fq_default_poly_add(d, d, b, ctx);
        FLINT_TEST(fq_default_poly_equal(a, d, ctx));

        fq_default_poly_set(a, b, ctx);
        fq_default_poly_scalar_addmul_fq_default(a, c, t, ctx);
        fq_default_poly_scalar_submul_fq_default(a, c, t, ctx);
        FLINT_TEST(fq_default_poly_equal(a, b, ctx));

        fq_default_poly_shift_left(a, b, 1, ctx);
        fq_default_poly_shift_right(a, a, 1, ctx);
        FLINT_TEST(fq_default_poly_equal(a, b, ctx));

        fq_default_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest(b, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest(c, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest(d, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest_not_zero(e, state, n_randint(state, 9) + 1, ctx);

        fq_default_poly_mul(a, a, e, ctx);
        fq_default_poly_mul(b, b, e, ctx);
        fq_default_poly_gcd(d, a, b, ctx);
        FLINT_TEST(fq_default_poly_divides(c, d, e, ctx));
        if (!fq_default_poly_is_zero(d, ctx))
        {
            FLINT_TEST(fq_default_poly_divides(c, a, d, ctx));
            fq_default_poly_mul(c, c, d, ctx);
            FLINT_TEST(fq_default_poly_equal(c, a, ctx));

            FLINT_TEST(fq_default_poly_divides(c, b, d, ctx));
            fq_default_poly_mul(c, c, d, ctx);
            FLINT_TEST(fq_default_poly_equal(c, b, ctx));
        }
        else
        {
            FLINT_TEST(fq_default_poly_is_zero(a, ctx));
            FLINT_TEST(fq_default_poly_is_zero(b, ctx));
        }

        fq_default_poly_xgcd(d, c, e, a, b, ctx);
        fq_default_poly_mul(c, c, a, ctx);
        fq_default_poly_mul(e, e, b, ctx);
        fq_default_poly_add(c, c, e, ctx);
        FLINT_TEST(fq_default_poly_equal(c, d, ctx));

        fq_default_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, ctx);
        fq_default_poly_randtest(c, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest(d, state, n_randint(state, 50), ctx);
        fq_default_poly_randtest(e, state, n_randint(state, 50), ctx);

        fq_default_poly_rem(e, a, b, ctx);
        fq_default_poly_divrem(c, d, a, b, ctx);
        FLINT_TEST(fq_default_poly_equal(d, e, ctx));
        fq_default_poly_mul(e, c, b, ctx);
        fq_default_poly_add(e, e, d, ctx);
        FLINT_TEST(fq_default_poly_equal(e, a, ctx));

        fq_default_clear(t, ctx);
        fq_default_poly_clear(a, ctx);
        fq_default_poly_clear(b, ctx);
        fq_default_poly_clear(c, ctx);
        fq_default_poly_clear(d, ctx);
        fq_default_poly_clear(e, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
