/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

/* check |res - ref| < 1.1 2^(e - prec), e the exponent of ref (computed with
   64 extra bits, which adds a negligible 2^-64 ulp) */
static int
check(const arf_t res, const arf_t ref, slong prec, arf_rnd_t rnd)
{
    arf_t d, bound;
    int ok;
    arf_init(d);
    arf_init(bound);
    arf_sub(d, res, ref, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_set_d(bound, (rnd == ARF_RND_ACCURATE) ? 0.51 : 1.0);
    arf_mul_2exp_fmpz(bound, bound, ARF_EXPREF(ref));
    arf_mul_2exp_si(bound, bound, -prec);
    ok = arf_cmpabs(d, bound) < 0 || arf_is_zero(d);
    arf_clear(d);
    arf_clear(bound);
    return ok;
}

#define CHECK(name) \
    if (!check(res, ref, prec, rnd)) \
        TEST_FUNCTION_FAIL("%s: prec = %wd\nx = %{arf}\ny = %{arf}\nres = %{arf}\nref = %{arf}\n", name, prec, x, y, res, ref);

TEST_FUNCTION_START(arf_approx, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, res, ref;
        slong prec, xbits, ybits;
        arf_rnd_t rnd = n_randint(state, 2) ? ARF_RND_FAST : ARF_RND_ACCURATE;

        arf_init(x); arf_init(y); arf_init(res); arf_init(ref);

        prec = 2 + n_randint(state, n_randint(state, 4) == 0 ? 300000 : 2000);
        xbits = 1 + n_randint(state, prec + 100);
        ybits = 1 + n_randint(state, prec + 100);
        if (n_randint(state, 2)) xbits = prec;
        if (n_randint(state, 2)) ybits = prec;

        arf_randtest_not_zero(x, state, xbits, 100);
        arf_randtest_not_zero(y, state, ybits, 100);
        if (n_randint(state, 4) == 0)
            arf_set_ui_2exp_si(y, 1 + n_randint(state, 100), n_randint(state, 200) - 100);

        arf_div(res, x, y, prec, rnd);
        arf_div(ref, x, y, prec + 64, ARF_RND_DOWN);
        CHECK("div")

        arf_ui_div(res, 1, y, prec, rnd);
        arf_ui_div(ref, 1, y, prec + 64, ARF_RND_DOWN);
        CHECK("inv")

        arf_abs(x, x);
        arf_sqrt(res, x, prec, rnd);
        arf_sqrt(ref, x, prec + 64, ARF_RND_DOWN);
        CHECK("sqrt")

        arf_rsqrt(res, x, prec, rnd);
        arf_rsqrt(ref, x, prec + 64, ARF_RND_DOWN);
        CHECK("rsqrt")

        /* aliased versions must agree with the unaliased ones */
        {
            arf_t t;
            arf_init(t);
            arf_div(res, x, y, prec, rnd);
            arf_set(t, x); arf_div(t, t, y, prec, rnd);
            if (!arf_equal(t, res)) TEST_FUNCTION_FAIL("div alias x: prec = %wd\n", prec);
            arf_set(t, y); arf_div(t, x, t, prec, rnd);
            if (!arf_equal(t, res)) TEST_FUNCTION_FAIL("div alias y: prec = %wd\n", prec);
            arf_ui_div(res, 1, y, prec, rnd); arf_set(t, y); arf_ui_div(t, 1, t, prec, rnd);
            if (!arf_equal(t, res)) TEST_FUNCTION_FAIL("inv alias: prec = %wd\n", prec);
            arf_sqrt(res, x, prec, rnd); arf_set(t, x); arf_sqrt(t, t, prec, rnd);
            if (!arf_equal(t, res)) TEST_FUNCTION_FAIL("sqrt alias: prec = %wd\n", prec);
            arf_rsqrt(res, x, prec, rnd); arf_set(t, x); arf_rsqrt(t, t, prec, rnd);
            if (!arf_equal(t, res)) TEST_FUNCTION_FAIL("rsqrt alias: prec = %wd\n", prec);
            arf_clear(t);
        }

        arf_clear(x); arf_clear(y); arf_clear(res); arf_clear(ref);
    }

    TEST_FUNCTION_END(state);
}
