/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acf.h"
#include "acb.h"

static void
_acb_set_acf(acb_t res, const acf_t x)
{
    arb_set_arf(acb_realref(res), acf_realref(x));
    arb_set_arf(acb_imagref(res), acf_imagref(x));
}

static int
_acf_is_zero(const acf_t x)
{
    return arf_is_zero(acf_realref(x)) && arf_is_zero(acf_imagref(x));
}

/* acf_approx_inv, acf_approx_div and acf_approx_sqrt against acb at a
   higher precision: the result within 2^(4 - prec) |ref| of the exact
   value (the kernels of arf promise a few ulps relative to the
   magnitude of the result), special operands (zero parts) included, and
   aliasing of the output with an operand */
TEST_FUNCTION_START(acf_approx_div, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        acf_t x, y, r;
        acb_t ax, ay, ar, az;
        mag_t m, t;
        slong prec, bits;
        int op, alias;
        arf_rnd_t rnd;

        acf_init(x);
        acf_init(y);
        acf_init(r);
        acb_init(ax);
        acb_init(ay);
        acb_init(ar);
        acb_init(az);
        mag_init(m);
        mag_init(t);

        prec = 2 + n_randint(state, n_randint(state, 4) ? 300 : 2000);
        bits = 1 + n_randint(state, n_randint(state, 4) ? 300 : 2000);
        acf_randtest(x, state, bits, 10);
        acf_randtest(y, state, bits, 10);
        if (n_randint(state, 8) == 0) arf_zero(acf_realref(x));
        if (n_randint(state, 8) == 0) arf_zero(acf_imagref(x));
        if (n_randint(state, 8) == 0) arf_zero(acf_realref(y));
        if (n_randint(state, 8) == 0) arf_zero(acf_imagref(y));
        if (n_randint(state, 8) == 0)
        {
            /* equal parts: cancellation in the quotient's numerator */
            arf_set(acf_imagref(y), acf_realref(y));
            if (n_randint(state, 2))
                arf_neg(acf_imagref(y), acf_imagref(y));
        }
        if (_acf_is_zero(y))
            arf_one(acf_realref(y));
        rnd = (arf_rnd_t) n_randint(state, 6);
        op = n_randint(state, 3);
        alias = n_randint(state, 3);

        _acb_set_acf(ax, x);
        _acb_set_acf(ay, y);

        if (op == 0)
        {
            if (_acf_is_zero(x))
                arf_one(acf_realref(x)), _acb_set_acf(ax, x);
            acb_inv(ar, ax, prec + 64);
            if (alias == 1)
            {
                acf_set(r, x);
                acf_approx_inv(r, r, prec, rnd);
            }
            else
                acf_approx_inv(r, x, prec, rnd);
        }
        else if (op == 1)
        {
            acb_div(ar, ax, ay, prec + 64);
            if (alias == 1)
            {
                acf_set(r, x);
                acf_approx_div(r, r, y, prec, rnd);
            }
            else if (alias == 2)
            {
                acf_set(r, y);
                acf_approx_div(r, x, r, prec, rnd);
            }
            else
                acf_approx_div(r, x, y, prec, rnd);
        }
        else
        {
            acb_sqrt(ar, ax, prec + 64);
            if (alias == 1)
            {
                acf_set(r, x);
                acf_approx_sqrt(r, r, prec, rnd);
            }
            else
                acf_approx_sqrt(r, x, prec, rnd);
        }

        if (!acb_is_finite(ar))
        {
            /* (division by zero and the like: nothing to check) */
        }
        else
        {
            _acb_set_acf(az, r);
            acb_sub(az, az, ar, prec + 64);
            acb_get_mag(m, az);
            acb_get_mag_lower(t, ar);
            mag_mul_2exp_si(t, t, 4 - prec);
            if (mag_cmp(m, t) > 0 && !(mag_is_zero(t) && mag_is_zero(m)))
            {
                flint_printf("FAIL: op %d, prec %wd, rnd %d, alias %d\n", op, prec, (int) rnd, alias);
                flint_printf("x = "); acf_printd(x, 20); flint_printf("\n");
                flint_printf("y = "); acf_printd(y, 20); flint_printf("\n");
                flint_printf("r = "); acf_printd(r, 20); flint_printf("\n");
                flint_printf("ref = "); acb_printn(ar, 20, 0); flint_printf("\n");
                flint_printf("|r - ref| = "); mag_printd(m, 5); flint_printf(", allowed "); mag_printd(t, 5); flint_printf("\n");
                flint_abort();
            }
        }

        acf_clear(x);
        acf_clear(y);
        acf_clear(r);
        acb_clear(ax);
        acb_clear(ay);
        acb_clear(ar);
        acb_clear(az);
        mag_clear(m);
        mag_clear(t);
    }

    TEST_FUNCTION_END(state);
}
