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

    flint_printf("jet_ellipsoid....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: sum of terms on border of ellipsoid must be less than bound */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        slong ord = n_randint(state, 3);
        slong nb = acb_theta_jet_nb(ord, g);
        acb_theta_eld_t E;
        acb_mat_t tau;
        acb_ptr z;
        acb_t term;
        arb_t u, abs, sum;
        slong nb_pts;
        slong* pts;
        slong* tups;
        slong k, j;

        acb_mat_init(tau, g, g);
        acb_theta_eld_init(E, g, g);
        z = _acb_vec_init(g);
        acb_init(term);
        arb_init(u);
        arb_init(abs);
        arb_init(sum);
        tups = flint_malloc(g * nb * sizeof(slong));

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g; k++)
        {
            acb_randtest_precise(&z[k], state, prec, bits);
        }

        /* Test: sum of terms on the border is less than u */
        acb_theta_jet_ellipsoid(E, u, z, tau, ord, prec);
        nb_pts = acb_theta_eld_nb_border(E);
        pts = flint_malloc(g * nb_pts * sizeof(slong));
        acb_theta_eld_border(pts, E);
        acb_theta_jet_tuples(tups, ord, g);

        for (j = 0; j < nb; j++)
        {
            arb_zero(sum);
            for (k = 0; k < nb_pts; k++)
            {
                acb_theta_naive_term(term, z, tau, tups + j * g, pts + k * g, 2 * prec);
                acb_abs(abs, term, 2 * prec);
                arb_add(sum, sum, abs, 2 * prec);
            }

            arb_sub(abs, sum, u, 2 * prec);

            if (!arb_is_negative(abs))
            {
                flint_printf("FAIL (ord = %wd, j = %wd)\n", ord, j);
                flint_printf("sum, bound:\n");
                arb_printd(sum, 10);
                flint_printf("\n");
                arb_printd(u, 10);
                flint_printf("\ntau:\n");
                acb_mat_printd(tau, 5);
                flint_printf("z:\n");
                _acb_vec_printd(z, g, 10);
                acb_theta_eld_print(E);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        acb_theta_eld_clear(E);
        _acb_vec_clear(z, g);
        acb_clear(term);
        arb_clear(u);
        arb_clear(abs);
        arb_clear(sum);
        flint_free(pts);
        flint_free(tups);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

