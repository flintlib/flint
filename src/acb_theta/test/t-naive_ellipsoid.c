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

    flint_printf("naive_ellipsoid....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: sum of terms on border of ellipsoid must be less than bound, and
       bounds are infinite on phony input */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        slong nbz = 1 + n_randint(state, 3);
        acb_theta_eld_t E, E2;
        acb_mat_t tau;
        acb_ptr c, z, new_z;
        arb_ptr u;
        acb_t term;
        arb_t abs, sum;
        slong nb_pts;
        slong* pts;
        slong k, j;

        acb_mat_init(tau, g, g);
        acb_theta_eld_init(E, g, g);
        acb_theta_eld_init(E2, g, g);
        z = _acb_vec_init(g * nbz);
        new_z = _acb_vec_init(g * nbz);
        c = _acb_vec_init(nbz);
        u = _arb_vec_init(nbz);
        acb_init(term);
        arb_init(abs);
        arb_init(sum);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g * nbz; k++)
        {
            acb_randtest_precise(&z[k], state, prec, bits);
        }

        /* Test: sum of terms on the border is less than u */
        acb_theta_naive_ellipsoid(E, new_z, c, u, z, nbz, tau, prec);
        nb_pts = acb_theta_eld_nb_border(E);
        pts = flint_malloc(g * nb_pts * sizeof(slong));
        acb_theta_eld_border(pts, E);

        for (j = 0; j < nbz; j++)
        {
            arb_zero(sum);
            for (k = 0; k < nb_pts; k++)
            {
                acb_theta_naive_term(term, new_z + j * g, tau, NULL, pts + k * g, 2 * prec);
                acb_abs(abs, term, 2 * prec);
                arb_add(sum, sum, abs, 2 * prec);
            }

            arb_sub(abs, sum, &u[j], 2 * prec);

            if (!arb_is_negative(abs))
            {
                flint_printf("FAIL\n");
                flint_printf("sum, bound:\n");
                arb_printd(sum, 10);
                flint_printf("\n");
                arb_printd(&u[j], 10);
                flint_printf("\ntau:\n");
                acb_mat_printd(tau, 5);
                flint_printf("new_z:\n");
                _acb_vec_printd(new_z + j * g, g, 10);
                acb_theta_eld_print(E);
                flint_abort();
            }
        }

        /* Test: indeterminate on phony tau */
        arb_randtest_positive(acb_imagref(acb_mat_entry(tau, 0, 0)), state, prec, bits);
        acb_neg(acb_mat_entry(tau, 0, 0), acb_mat_entry(tau, 0, 0));
        acb_theta_naive_ellipsoid(E2, u, z, tau, ord, prec);
        if (acb_is_finite(c) && arb_is_finite(u))
        {
            flint_printf("FAIL (not infinite)\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_theta_eld_clear(E);
        acb_theta_eld_clear(E2);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(new_z, g * nbz);
        _acb_vec_clear(c, nbz);
        _arb_vec_clear(u, nbz);
        acb_clear(term);
        arb_clear(abs);
        arb_clear(sum);
        flint_free(pts);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
