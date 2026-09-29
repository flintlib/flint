/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "gr.h"

/* random rational with numerator/denominator up to bits */
static void
_random_fmpq(fmpq_t q, flint_rand_t state, slong bits)
{
    fmpz_randtest(fmpq_numref(q), state, bits);
    do {
        fmpz_randtest_not_zero(fmpq_denref(q), state, bits);
    } while (fmpz_is_zero(fmpq_denref(q)));
    fmpz_abs(fmpq_denref(q), fmpq_denref(q));

    /* sometimes make it a decimal fraction or an integer */
    switch (n_randint(state, 4))
    {
        case 0:
            fmpz_ui_pow_ui(fmpq_denref(q), 10, n_randint(state, 40));
            break;
        case 1:
            fmpz_one(fmpq_denref(q));
            break;
        case 2:
            fmpz_ui_pow_ui(fmpq_denref(q), 2, n_randint(state, 40));
            break;
        default:
            break;
    }

    fmpq_canonicalise(q);
}

TEST_FUNCTION_START(decfloat_round, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y, z;
        fmpq_t q;
        slong prec;
        int rnd, s1, s2;

        gr_ctx_init_decfloat_randtest(ctx, state, 40);
        /* no exponent limits for this test */
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(z, ctx);
        fmpq_init(q);

        prec = 1 + n_randint(state, 60);
        if (n_randint(state, 10) == 0)
            prec = DECIMAL_PREC_EXACT;
        rnd = n_randint(state, DECIMAL_RND_NUM);

        /* rounding of rationals */
        _random_fmpq(q, state, 1 + n_randint(state, 200));

        s1 = decfloat_set_round_fmpq(x, q, prec, rnd, ctx);
        s2 = decfloat_set_round_fmpq_reference(y, q, prec, rnd, ctx);

        if (s1 != s2 || (s1 == GR_SUCCESS && decfloat_equal(x, y, ctx) != T_TRUE))
        {
            flint_printf("FAIL: set_round_fmpq\n");
            gr_ctx_println(ctx);
            flint_printf("prec = %wd, rnd = %d\n", prec, rnd);
            flint_printf("q = %{fmpq}\n", q);
            flint_printf("x = %{gr} (%d)\n", x, ctx, s1);
            flint_printf("y = %{gr} (%d)\n", y, ctx, s2);
            flint_abort();
        }

        /* re-rounding an already rounded value to a lower precision */
        if (s1 == GR_SUCCESS && prec != DECIMAL_PREC_EXACT)
        {
            slong prec2 = 1 + n_randint(state, prec);
            int rnd2 = n_randint(state, DECIMAL_RND_NUM);
            fmpq_t q2;

            fmpq_init(q2);

            if (decfloat_get_fmpq(q2, x, ctx) == GR_SUCCESS)
            {
                if (n_randint(state, 2))
                {
                    s1 = decfloat_set_round(z, x, prec2, rnd2, ctx);
                }
                else
                {
                    GR_MUST_SUCCEED(decfloat_set_round(z, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
                    s1 = decfloat_set_round(z, z, prec2, rnd2, ctx);
                }

                s2 = decfloat_set_round_fmpq_reference(y, q2, prec2, rnd2, ctx);

                if (s1 != s2 || (s1 == GR_SUCCESS && decfloat_equal(z, y, ctx) != T_TRUE))
                {
                    flint_printf("FAIL: set_round\n");
                    gr_ctx_println(ctx);
                    flint_printf("prec2 = %wd, rnd2 = %d\n", prec2, rnd2);
                    flint_printf("x = %{gr}\n", x, ctx);
                    flint_printf("z = %{gr} (%d)\n", z, ctx, s1);
                    flint_printf("y = %{gr} (%d)\n", y, ctx, s2);
                    flint_abort();
                }
            }

            fmpq_clear(q2);
        }

        /* significant digit count is at most prec, and the representation is canonical */
        if (s1 == GR_SUCCESS && !DECFLOAT_IS_SPECIAL(x))
        {
            slong n = FLINT_ABS(x->m.size);

            if (x->m.d[0] == 0 || x->m.d[n - 1] == 0)
            {
                flint_printf("FAIL: non-canonical mantissa\n");
                flint_abort();
            }

            if (prec != DECIMAL_PREC_EXACT && decfloat_digits(x, ctx) > prec)
            {
                flint_printf("FAIL: too many significant digits\n");
                flint_printf("x = %{gr}, prec = %wd, sd = %wd\n", x, ctx, prec, decfloat_digits(x, ctx));
                flint_abort();
            }
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(z, ctx);
        fmpq_clear(q);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
