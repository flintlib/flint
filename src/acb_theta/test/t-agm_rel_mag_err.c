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

    flint_printf("agm_rel_mag_err....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: recover expected values */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, 8);
        slong prec = 100 + n_randint(state, 1000);
        slong hprec = prec + 100;
        slong bits = n_randint(state, 5);
        acb_ptr a;
        arb_ptr dist;
        arf_t m, eps, t_m, t_eps;
        arb_t x;
        slong k;

        a = _acb_vec_init(n);
        dist = _arb_vec_init(n);
        arf_init(m);
        arf_init(eps);
        arf_init(t_m);
        arf_init(t_eps);
        arb_init(x);

        /* Generate m, eps, dist */
        arf_randtest(m, state, prec, bits);
        if (arf_cmp_si(m, 0) < 0)
        {
            arf_neg(m, m);
        }
        arf_one(eps);
        arf_mul_2exp_si(eps, eps, -prec);
        for (k = 0; k < n; k++)
        {
            arb_randtest_positive(&dist[k], state, prec, bits);
        }

        /* Generate values */
        for (k = 0; k < n; k++)
        {
            if (k == 0)
            {
                acb_one(&a[k]);
            }
            else
            {
                acb_urandom(&a[k], state, hprec);
            }
            arb_neg(x, &dist[k]);
            arb_exp(x, x, prec);
            arb_mul_arf(x, x, m, hprec);
            acb_mul_arb(&a[k], &a[k], x, hprec);

            arb_neg(x, &dist[k]);
            arb_exp(x, x, prec);
            arb_mul_arf(x, x, eps, hprec);
            acb_add_error_arb(&a[k], x);
        }

        acb_theta_agm_rel_mag_err(t_m, t_eps, a, dist, n, prec);

        if (arf_cmp(t_m, m) < 0 || arf_cmp(t_eps, eps) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("n = %wd, distances:\n", n);
            _arb_vec_printd(dist, n, 5);
            flint_printf("values:\n");
            _acb_vec_printd(a, n, 5);
            flint_printf("m, eps, t_m, t_eps: ");
            arf_printd(m, 5);
            flint_printf(", ");
            arf_printd(eps, 5);
            flint_printf(", ");
            arf_printd(t_m, 5);
            flint_printf(", ");
            arf_printd(t_eps, 5);
            flint_printf("\n");
            flint_abort();
        }

        _acb_vec_clear(a, n);
        _arb_vec_clear(dist, n);
        arf_clear(m);
        arf_clear(eps);
        arf_clear(t_m);
        arf_clear(t_eps);
        arb_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

