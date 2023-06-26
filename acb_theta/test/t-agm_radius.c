/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: correct radius after a few good steps */
    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 500 + n_randint(state, 1000);
        slong mag_bits = 1 + n_randint(state, 4);
        slong nb_good = 1 + n_randint(state, 10);
        slong quo = 2 + n_randint(state, 10);
        slong lowprec = ACB_THETA_AGM_LOWPREC;
        slong n = 1 << g;
        acb_ptr a;
        acb_ptr b;
        acb_t x;
        arb_t abs;
        arf_t rad;
        arf_t test;
        arf_struct *mi;
        arf_struct *Mi;
        slong k;

        a = _acb_vec_init(n);
        b = _acb_vec_init(n);
        acb_init(x);
        arb_init(abs);
        arf_init(rad);
        arf_init(test);
        mi = flint_malloc(nb_good * sizeof(arf_struct));
        Mi = flint_malloc(nb_good * sizeof(arf_struct));
        for (k = 0; k < nb_good; k++)
        {
            arf_init(&mi[k]);
            arf_init(&Mi[k]);
        }

        /* Generate starting vector */
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, -1);

        acb_one(x);
        for (k = 0; k < n; k++)
            acb_randtest_disk(&a[k], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a, a, n, x, prec);

        /* Get mi, Mi */
        _acb_vec_set(b, a, n);
        for (k = 0; k < nb_good; k++)
        {
            acb_theta_agm_max_abs(abs, b, n, lowprec);
            arb_get_ubound_arf(&Mi[k], abs, lowprec);
            acb_theta_agm_min_abs(abs, b, n, lowprec);
            arb_get_lbound_arf(&mi[k], abs, lowprec);
            acb_theta_agm_step_good(b, b, g, prec);
        }

        /* Choose epsilon, compute rad */
        acb_theta_agm_abs_dist(abs, a, n, lowprec, prec);
        arb_div_si(abs, abs, quo, lowprec);
        arb_get_ubound_arf(rad, abs, lowprec);

        acb_theta_agm_radius(test, mi, Mi, rad, nb_good, prec);

        /* Generate perturbed data, compute steps, check distance */
        for (k = 0; k < n; k++)
            acb_randtest_disk(&b[k], x, test, state, prec);
        _acb_vec_add(b, b, a, n, prec);
        for (k = 0; k < nb_good; k++)
            acb_theta_agm_step_good(b, b, g, prec);
        acb_theta_agm_abs_dist(abs, b, n, lowprec, prec);
        arb_get_lbound_arf(test, abs, lowprec);

        if (arf_cmp(test, rad) > 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        _acb_vec_clear(a, n);
        _acb_vec_clear(b, n);
        acb_clear(x);
        arb_clear(abs);
        arf_clear(rad);
        arf_clear(test);
        for (k = 0; k < nb_good; k++)
        {
            arf_clear(&mi[k]);
            arf_clear(&Mi[k]);
        }
        flint_free(mi);
        flint_free(Mi);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
