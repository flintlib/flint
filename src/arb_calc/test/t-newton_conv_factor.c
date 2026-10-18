/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "arb_calc.h"

/* f(x) = x^2 - c, with Taylor coefficients f, f', f''/2 */
static int
quadratic(arb_ptr out, const arb_t x, void * param, slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_mul(t, x, x, prec);
    if (order > 1)
        arb_mul_2exp_si(out + 1, x, 1);
    if (order > 2)
        arb_one(out + 2);
    if (order > 3)
        _arb_vec_zero(out + 3, order - 3);
    arb_sub(out, t, (arb_srcptr) param, prec);
    arb_clear(t);
    return 0;
}

TEST_FUNCTION_START(arb_calc_newton_conv_factor, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t c, s, I, x, y, t;
        arf_t C, lo;
        slong prec;
        ulong cn;
        int r;

        arb_init(c); arb_init(s); arb_init(I); arb_init(x); arb_init(y);
        arb_init(t);
        arf_init(C); arf_init(lo);

        prec = 30 + n_randint(state, 200);

        /* the root s = sqrt(c) for c = cn / 4 not a square */
        do {
            cn = 2 + n_randint(state, 1000);
        } while (n_is_square(cn));
        arb_set_ui(c, cn);
        arb_mul_2exp_si(c, c, -2);
        arb_sqrt(s, c, 2 * prec + 100);

        /* conv region I = [s +/- 2^-k], start x = [s + d +/- 2^-k / 4]
           with |d| just below 2^-k / 4 */
        {
            slong k = 2 + n_randint(state, 10);
            arb_set(I, s);
            mag_set_ui_2exp_si(arb_radref(I), 1, -k);
            arb_set_round(I, I, 64);

            arb_set_ui(t, 255 - n_randint(state, 8));
            arb_mul_2exp_si(t, t, -k - 10);
            if (n_randint(state, 2))
                arb_neg(t, t);
            arb_add(x, s, t, prec);
            arb_get_mid_arb(x, x);
            mag_set_ui_2exp_si(arb_radref(x), 1, -k - 2);
        }

        if (!arb_contains(x, s) || !arb_contains(I, x))
            TEST_FUNCTION_FAIL("bad test input\n");

        /* C must bound sup |f''| / (2 inf |f'|) = 1 / (2 inf I) */
        arb_calc_newton_conv_factor(C, quadratic, c, I, prec);
        arb_get_lbound_arf(lo, I, ARF_PREC_EXACT);
        arf_mul_2exp_si(lo, lo, 1);
        arf_ui_div(lo, 1, lo, prec, ARF_RND_DOWN);
        if (arf_cmp(C, lo) < 0)
            TEST_FUNCTION_FAIL("conv_factor\nI = %{arb}\nC = %{arf}\nlower bound = %{arf}\n", I, C, lo);

        r = arb_calc_newton_step(y, quadratic, c, x, I, C, prec);
        if (!arb_contains(y, s))
            TEST_FUNCTION_FAIL("newton_step (status %d)\nx = %{arb}\ny = %{arb}\ns = %{arb}\n", r, x, y, s);

        r = arb_calc_refine_root_newton(y, quadratic, c, x, I, C, 10, 2 * prec);
        if (!arb_contains(y, s))
            TEST_FUNCTION_FAIL("refine_root_newton (status %d)\nx = %{arb}\ny = %{arb}\ns = %{arb}\n", r, x, y, s);

        arb_clear(c); arb_clear(s); arb_clear(I); arb_clear(x); arb_clear(y);
        arb_clear(t);
        arf_clear(C); arf_clear(lo);
    }

    /* the example from the issue: f(x) = x^2 - 2 on [1.4, 1.5] */
    {
        arb_t c, s, I, x, y;
        arf_t C;

        arb_init(c); arb_init(s); arb_init(I); arb_init(x); arb_init(y);
        arf_init(C);

        arb_set_ui(c, 2);
        arb_sqrt_ui(s, 2, 200);
        arb_set_str(I, "[1.45 +/- 0.05]", 53);
        arb_calc_newton_conv_factor(C, quadratic, c, I, 53);
        if (arf_cmp_d(C, 0.357) < 0)
            TEST_FUNCTION_FAIL("conv_factor on [1.4, 1.5]: C = %{arf}\n", C);
        arb_set_str(x, "[1.415 +/- 0.001]", 200);
        arb_calc_newton_step(y, quadratic, c, x, I, C, 200);
        if (!arb_contains(y, s))
            TEST_FUNCTION_FAIL("newton_step from [1.415 +/- 0.001]: y = %{arb}\n", y);

        arb_clear(c); arb_clear(s); arb_clear(I); arb_clear(x); arb_clear(y);
        arf_clear(C);
    }

    TEST_FUNCTION_END(state);
}
