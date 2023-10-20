/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_sqrt....");
    fflush(stdout);

    flint_randinit(state);

    /* Test:
       - if nonzero, value of square root should agree; precision remains high
       - if contains zero or random values, should contain both square roots */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_t rt;
        acb_t x;
        arb_t err;
        acb_t rt_low;
        acb_t test;

        slong prec = 100 + n_randint(state, 1000);
        slong mag_bits = n_randint(state, 4);
        slong lowprec = 10 + n_randint(state, 10);

        acb_init(rt);
        acb_init(x);
        arb_init(err);
        acb_init(rt_low);
        acb_init(test);

        acb_randtest_precise(rt, state, prec, mag_bits);
        if (iter % 10 == 0)
        {
            acb_union(rt, rt, x, prec);
        }

        acb_sqr(x, rt, prec);
        arb_one(err);
        arb_mul_2exp_si(err, err, -lowprec);
        arb_add_si(err, err, 1, lowprec);
        acb_mul_arb(rt_low, rt, err, lowprec);

        if (iter % 10 == 1)
        {
            acb_randtest(rt_low, state, prec, mag_bits);
        }

        acb_theta_agm_sqrt(test, x, rt_low, 1, prec);

        if (!acb_contains(test, rt))
        {
            flint_printf("FAIL (value)\n");
            fflush(stdout);
            flint_abort();
        }

        if (acb_contains(rt_low, rt) && !acb_is_finite(test))
        {
            flint_printf("FAIL (infinite)\n");
            fflush(stdout);
            flint_abort();
        }

        acb_get_mid(x, test);
        acb_sub(test, test, x, prec);
        acb_abs(err, test, prec);
        arb_mul_2exp_si(err, err, prec - n_pow(2, mag_bits) - 10);
        arb_add_si(err, err, -1, prec);

        if (!arb_contains_zero(rt) && !arb_is_negative(err))
        {
            flint_printf("FAIL (precision)\n");
            flint_printf("prec = %wd, mag_bits = %wd, difference:\n", prec, mag_bits);
            acb_printd(test, 10);
            flint_printf("\n");
            flint_printf("rt_low:\n");
            acb_printd(rt_low, 10);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        acb_clear(rt);
        acb_clear(x);
        arb_clear(err);
        acb_clear(rt_low);
        acb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
