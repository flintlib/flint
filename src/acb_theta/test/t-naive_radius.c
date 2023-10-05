/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Evaluate upper bound on the tail */
static void
acb_theta_naive_tail(arb_t res, const arf_t R2, const arb_mat_t C, slong ord)
{
    slong g = arb_mat_nrows(C);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t t, Rm;
    slong k;

    arb_init(t);
    arb_init(Rm);

    /* Ensure assumptions R2\geq 4, R2\geq 2*ord are satisfied */
    arb_set_arf(Rm, R2);
    arb_set_si(t, FLINT_MAX(4, 2 * ord));
    arb_max(Rm, Rm, t, lp);

    /* Evaluate 2^(2*g+2) R^(g - 1 + ord) exp(-R^2) \prod(1 + gamma_i^{-1}) */
    arb_one(res);
    arb_mul_2exp_si(res, res, 2 * g + 2);

    arb_sqrt(t, Rm, lp);
    arb_pow_ui(t, t, g - 1 + ord, lp);
    arb_mul(res, res, t, lp);

    arb_neg(t, Rm);
    arb_exp(t, t, lp);
    arb_mul(res, res, t, lp);

    for (k = 0; k < g; k++)
    {
        arb_inv(t, arb_mat_entry(C, k, k), lp);
        arb_add_si(t, t, 1, lp);
        arb_mul(res, res, t, lp);
    }

    arb_clear(t);
    arb_clear(Rm);
}

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: tail is small */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong ord = n_randint(state, 10);
        slong prec = ACB_THETA_LOW_PREC;
        slong bits = n_randint(state, 5);
        slong exp = 10 + n_randint(state, 100);
        arb_mat_t Y;
        arf_t R2, eps, t;
        arb_t bound;

        arb_mat_init(Y, g, g);
        arf_init(R2);
        arf_init(eps);
        arf_init(t);
        arb_init(bound);

        arb_mat_randtest_cho(Y, state, prec, bits);
        arb_mat_transpose(Y, Y);

        acb_theta_naive_radius(R2, eps, Y, ord, exp);
        acb_theta_naive_tail(bound, R2, Y, ord);
        arb_get_lbound_arf(t, bound, prec);

        if (arf_cmp(t, eps) > 0)
        {
            flint_printf("FAIL\n");
            arb_mat_printd(Y, 5);
            flint_printf("exp = %wd, ord = %wd, eps, R2:\n", exp, ord);
            arf_printd(eps, 10);
            flint_printf("\n");
            arf_printd(R2, 10);
            flint_printf("\nbound:\n");
            arb_printd(bound, 10);
            flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(Y);
        arf_clear(R2);
        arf_clear(eps);
        arf_clear(t);
        arb_clear(bound);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
