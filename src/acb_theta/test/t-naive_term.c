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

    flint_printf("naive_term....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with genus 1 */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1;
        slong prec = 100 + n_randint(state, 200);
        slong bits = n_randint(state, 5);
        slong n = n_randint(state, 100);
        slong k = n_randint(state, 10);
        acb_mat_t tau;
        acb_t z;
        acb_t x, t;

        acb_mat_init(tau, g, g);
        acb_init(z);
        acb_init(x);
        acb_init(t);

        acb_siegel_randtest(tau, state, prec, bits);
        acb_randtest_precise(z, state, prec, bits);

        acb_theta_naive_term(x, z, tau, &k, &n, prec);
        acb_mul_si(t, acb_mat_entry(tau, 0, 0), n, prec);
        acb_addmul_si(t, z, 2, prec);
        acb_mul_si(t, t, n, prec);
        acb_exp_pi_i(t, t, prec);

        fmpz_pow_ui(n, n, k);
        acb_mul(t, t, n, prec);

        if (!acb_overlaps(x, t))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_clear(z);
        acb_clear(x);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
