/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

/* old acb_log code for comparison */
static int
close_to_one(const acb_t z)
{
    mp_limb_t top;

    if (arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z))) > -3)
        return 0;

    if (ARF_EXP(arb_midref(acb_realref(z))) == 0)
    {
        ARF_GET_TOP_LIMB(top, arb_midref(acb_realref(z)));

        return (top >> (FLINT_BITS - 4)) == 15;
    }
    else if (ARF_EXP(arb_midref(acb_realref(z))) == 1)
    {
        ARF_GET_TOP_LIMB(top, arb_midref(acb_realref(z)));

        return (top >> (FLINT_BITS - 4)) == 8;
    }

    return 0;
}

void
acb_log_old(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)

    if (arb_is_zero(b))
    {
        if (arb_is_positive(a))
        {
            arb_log(acb_realref(r), a, prec);
            arb_zero(acb_imagref(r));
        }
        else if (arb_is_negative(a))
        {
            arb_neg(acb_realref(r), a);
            arb_log(acb_realref(r), acb_realref(r), prec);
            arb_const_pi(acb_imagref(r), prec);
        }
        else
        {
            acb_indeterminate(r);
        }
    }
    else if (arb_is_zero(a))
    {
        if (arb_is_positive(b))
        {
            arb_log(acb_realref(r), b, prec);
            arb_const_pi(acb_imagref(r), prec);
            arb_mul_2exp_si(acb_imagref(r), acb_imagref(r), -1);
        }
        else if (arb_is_negative(b))
        {
            arb_neg(acb_realref(r), b);
            arb_log(acb_realref(r), acb_realref(r), prec);
            arb_const_pi(acb_imagref(r), prec);
            arb_mul_2exp_si(acb_imagref(r), acb_imagref(r), -1);
            arb_neg(acb_imagref(r), acb_imagref(r));
        }
        else
        {
            acb_indeterminate(r);
        }
    }
    else
    {
        arb_t t, u;

        arb_init(t);
        arb_init(u);

        if (close_to_one(z))
        {
            arb_sub_ui(u, a, 1, prec + 8);
            arb_mul(t, u, u, prec + 8);
            arb_addmul(t, b, b, prec + 8);
            arb_mul_2exp_si(u, u, 1);
            arb_add(t, t, u, prec + 8);

            arb_log1p(t, t, prec);
            arb_mul_2exp_si(t, t, -1);
        }
        else
        {
            arb_mul(t, a, a, prec + 8);
            arb_addmul(t, b, b, prec + 8);

            if (arb_contains_zero(t) || arf_sgn(arb_midref(t)) < 0)
                arb_zero_pm_inf(t);
            else
                arb_log(t, t, prec);

            arb_mul_2exp_si(t, t, -1);
        }

        acb_arg(u, z, prec);

        arb_swap(acb_realref(r), t);
        arb_swap(acb_imagref(r), u);

        arb_clear(t);
        arb_clear(u);
    }

    if (!acb_is_finite(r))
        acb_indeterminate(r);

#undef a
#undef b
}

TEST_FUNCTION_START(acb_log, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, a, b, d;
        slong prec1, prec2, acc1, acc2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(a);
        acb_init(b);
        acb_init(d);

        acb_randtest_special(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest_special(a, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest_special(b, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            acb_log(a, x, prec1);
        }
        else
        {
            acb_set(a, x);  /* test aliasing */
            acb_log(a, a, prec1);
        }

        if (n_randint(state, 2))
        {
            acb_log(b, x, prec2);
        }
        else
        {
            acb_set(b, x);  /* test aliasing */
            acb_log(b, b, prec2);
        }

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_contains_zero(x) && acb_is_finite(x) && !acb_is_finite(a))
        {
            flint_printf("FAIL: not finite\n\n");
            flint_printf("prec1 = %wd\n\n", prec1);
            flint_printf("x = "); acb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");

            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");

            flint_abort();
        }

        /* check exp(log(x)) = x */
        acb_exp(b, b, prec1);

        if (!acb_contains(b, x))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_log_old(d, x, prec1);

        /* check consistency */
        if (!acb_overlaps(a, d))
        {
            flint_printf("FAIL: overlap 2\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* compare accuracy with log_old */
        acc1 = arb_rel_accuracy_bits(acb_realref(a));
        acc2 = arb_rel_accuracy_bits(acb_realref(d));

        if (acc2 > 0 && acc1 < acc2 - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
            flint_printf("x = "); acb_printd(x, 50); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 50); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 50); flint_printf("\n\n");

            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(acb_imagref(a));
        acc2 = arb_rel_accuracy_bits(acb_imagref(d));

        if (acc2 > 0 && acc1 < acc2 - 2)
        {
            flint_printf("FAIL: accuracy 2\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
            flint_printf("x = "); acb_printd(x, 50); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 50); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 50); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(a);
        acb_clear(b);
        acb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
