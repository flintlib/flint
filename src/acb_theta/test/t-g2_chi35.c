/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
chi35_lead(acb_t r, const acb_mat_t tau, slong prec)
{
    acb_t q1, q2, q3, t;

    acb_init(q1);
    acb_init(q2);
    acb_init(q3);
    acb_init(t);

    acb_exp_pi_i(q1, acb_mat_entry(tau, 0, 0), prec);
    acb_exp_pi_i(q2, acb_mat_entry(tau, 0, 1), prec);
    acb_exp_pi_i(q3, acb_mat_entry(tau, 1, 1), prec);

    acb_mul_(r, q1, q3, prec);
    acb_sqr(r, r, prec);
    acb_sub(t, q1, q3, prec);
    acb_mul(r, r, t, prec);
    acb_inv(t, q2, prec);
    acb_sub(t, q2, t, prec);
    acb_mul(r, r, t, prec);

    acb_clear(q1);
    acb_clear(q2);
    acb_clear(q3);
    acb_clear(t);
}

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("g2_chi35....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: transforms like a modular form */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n2 = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 10);
        fmpz_mat_t mat;
        acb_mat_t tau;
        acb_ptr th, z;
        slong k;
        acb_t r, s;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        acb_mat_init(tau, g, g);
        th = _acb_vec_init(n2);
        z = _acb_vec_init(g);
        acb_init(r);
        acb_init(s);

        sp2gz_randtest(mat, state, mag_bits);
        for (k = 0; k < n2; k++)
        {
            acb_randtest_precise(&th[k], state, prec, mag_bits);
        }

        acb_theta_g2_chi35(r, th, prec);
        acb_theta_transform_proj(th, mat, th, 0, prec);
        acb_theta_g2_chi35(s, th, prec);
        acb_mul_powi(s, s, acb_theta_transform_kappa(mat));

        if (!acb_overlaps(r, s))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(s, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_siegel_randtest_nice(tau, state, prec);
        acb_mul_2exp_si(tau, tau, n_randint(state, 10));
        acb_theta_naive_all(th, z, 1, tau, prec);
        acb_theta_g2_chi35(r, th, prec);
        chi35_lead(s, tau, prec);

        flint_printf("INFO: value/lead for matrix\n");
        acb_mat_printd(tau, 5);
        acb_printd(r, 10);
        flint_printf("\n");
        acb_printd(s, 10);
        flint_printf("\n");

        fmpz_mat_clear(mat);
        acb_mat_clear(tau);
        _acb_vec_clear(th, n2);
        _acb_vec_clear(z, g);
        acb_clear(r);
        acb_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
