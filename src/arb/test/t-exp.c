/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "fmpq.h"
#include "arb.h"

/* check accuracy compared to reference algorithm */
void arb_exp_simple(arb_t res, const arb_t x, slong prec)
{
    mag_t t, u;

    mag_init_set(t, arb_radref(x));
    mag_init(u);

    arf_set(arb_midref(res), arb_midref(x));
    mag_zero(arb_radref(res));
    arb_exp(res, x, prec);

    mag_expm1(t, t);
    arb_get_mag(u, res);
    mag_addmul(arb_radref(res), t, u);

    mag_clear(t);
    mag_clear(u);
}

TEST_FUNCTION_START(arb_exp, state)
{
    slong iter;

    /* check large arguments + compare with exp_simple */
    for (iter = 0; iter < 100000 *0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec0, prec1, prec2, acc1, acc2;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);

        arb_exp(b, a, prec1);
        arb_exp(c, a, prec2);
        arb_exp_simple(d, a, prec1);

        if (!arb_overlaps(b, c) || !arb_overlaps(b, d) || !arb_overlaps(c, d))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        /* compare accuracy with exp_simple */
        acc1 = arb_rel_accuracy_bits(b);
        acc2 = arb_rel_accuracy_bits(d);

        if (acc2 > 0 && acc1 < acc2 - 1)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
            flint_printf("a = "); arb_printd(a, 50); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 50); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);

        /* check exp(a)*exp(b) = exp(a+b) */
        arb_exp(c, a, prec1);
        arb_exp(d, b, prec1);
        arb_mul(c, c, d, prec1);

        arb_add(d, a, b, prec1);
        arb_exp(d, d, prec1);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    /* test union */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d, e;
        slong prec0, prec1, prec2, prec3, prec4;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);
        prec3 = 2 + n_randint(state, prec0);
        prec4 = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 200);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 200);
        arb_randtest_special(c, state, 1 + n_randint(state, prec0), 200);
        arb_randtest_special(d, state, 1 + n_randint(state, prec0), 200);
        arb_randtest_special(e, state, 1 + n_randint(state, prec0), 200);

        arb_exp(c, a, prec1);
        arb_exp(d, b, prec2);
        arb_union(e, a, b, prec3);
        arb_exp(e, e, prec4);

        if (!arb_overlaps(e, c) || !arb_overlaps(e, d))
        {
            flint_printf("FAIL: union\n\n");
            flint_printf("a = "); arb_printn(a, 1000, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 1000, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 1000, 0); flint_printf("\n\n");
            flint_printf("d = "); arb_printn(d, 1000, 0); flint_printf("\n\n");
            flint_printf("e = "); arb_printn(e, 1000, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(e);
    }

    /* comparison with mpfr */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        slong prec0, prec;

        prec0 = 400;
        if (iter % 100 == 0)
            prec0 = 10000;

        prec = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, prec0), 4);
        arb_randtest(b, state, 1 + n_randint(state, prec0), 4);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, prec0));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_exp(t, t, MPFR_RNDN);

        arb_exp(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("iter = %wd, prec = %wd\n\n", iter, prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_exp(a, a, prec);

        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("iter = %wd, prec = %wd\n\n", iter, prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    TEST_FUNCTION_END(state);
}
