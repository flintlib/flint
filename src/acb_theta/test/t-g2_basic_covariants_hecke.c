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

    flint_printf("g2_basic_covariants_hecke....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: find eigenvalue (1+p^2)(1+p^3) for E4 = -(Co20^2 - 3*Co40)/20 at p prime */
    for (iter = 0; iter < 1 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong nb_cov = ACB_THETA_G2_BASIC_NB;
        acb_mat_t tau;
        acb_ptr cov;
        acb_t r, s, t, u, v;
        slong prec = 100;
        slong primes[] = {2,3,5}; /*,3,5,7,11,13,17};*/
        slong nprimes = 3;
        slong k, p, l;

        acb_mat_init(tau, g, g);
        acb_init(u);
        acb_init(v);
        acb_init(r);
        acb_init(s);
        acb_init(t);

        /* Generate matrix with nice imaginary part */
        arb_urandom(acb_realref(acb_mat_entry(tau, 0, 0)), state, prec);
        arb_set_si(acb_imagref(acb_mat_entry(tau, 0, 0)), 1);
        arb_set_si(acb_imagref(acb_mat_entry(tau, 1, 1)), 1);
        arb_urandom(acb_realref(acb_mat_entry(tau, 0, 1)), state, prec);
        arb_urandom(acb_imagref(acb_mat_entry(tau, 0, 1)), state, prec);
        acb_mul_2exp_si(acb_mat_entry(tau, 0, 1), acb_mat_entry(tau, 0, 1), -2);
        acb_set(acb_mat_entry(tau, 1, 0), acb_mat_entry(tau, 0, 1));

        for (k = 0; k < nprimes; k++)
        {
            p = primes[k];
            /* flint_printf("\n\n\n*** Start p = %wd ***\n\n", p); */

            cov = _acb_vec_init(nb_cov * (acb_theta_g2_hecke_nb(p) + 1));
            acb_theta_g2_basic_covariants_hecke(cov, tau, p, prec);

            /* Get Co20 - 3*Co40 at tau */
            acb_mul_si(u, &cov[8], -3, prec);
            acb_sqr(v, &cov[1], prec);
            acb_add(s, u, v, prec);

            /* Get sum of Co20 - 3*Co40 at images */
            acb_zero(r);
            for (l = 0; l < acb_theta_g2_hecke_nb(p); l++)
            {
                acb_mul_si(u, &cov[(l + 1) * nb_cov + 8], -3, prec);
                acb_sqr(v, &cov[(l + 1) * nb_cov + 1], prec);
                acb_add(u, u, v, prec);
                acb_add(r, r, u, prec);
            }

            acb_div(r, r, s, prec);
            acb_set_si(s, n_pow(p, 5));
            acb_mul(r, r, s, prec);

            /* Get expected eigenvalue */
            if (n_is_prime(p))
            {
                acb_set_si(t, (1 + n_pow(p, 2)) * (1 + n_pow(p, 3)));
            }
            else
            {
                p = n_sqrt(p);
                acb_set_si(t, n_pow(p, 4) - n_pow(p, 2) + n_pow(p, 6)
                    + n_pow(p, 7) + p + n_pow(p, 2));
            }

            if (!acb_overlaps(r, t))
            {
                flint_printf("FAIL (p = %wd)\n", p);
                acb_printd(r, 5);
                flint_printf("\n");
                acb_printd(t, 5);
                flint_printf("\n");
                flint_abort();
            }

            _acb_vec_clear(cov, nb_cov * (acb_theta_g2_hecke_nb(p) + 1));
        }

        acb_mat_clear(tau);
        acb_clear(u);
        acb_clear(v);
        acb_clear(r);
        acb_clear(s);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

