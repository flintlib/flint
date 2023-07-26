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

    flint_printf("transform_sqr....");
    fflush(stdout);

    /* Test: compatible with transformation in genus 1 */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1;
        slong n2 = 1 << (2 * g);
        slong prec = 100;
        slong bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z;
        fmpz_mat_t mat;
        acb_ptr th, test;
        slong k2;
        slong k;

        acb_mat_init(tau, g, g);
        fmpz_mat_init(mat, 2 * g, 2 * g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(n2);
        test = _acb_vec_init(n2);

        acb_siegel_randtest_nice(tau, state, prec);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        sp2gz_randtest(mat, state, bits);
        k2 = acb_theta_transform_k2(mat);

        acb_theta_ql_all_sqr(th, z, tau, prec);
        acb_theta_transform_sqr(th, th, z, tau, mat, k2, prec);

        acb_siegel_transform_ext(z, tau, mat, z, tau, prec);
        acb_modular_theta(&test[3], &test[2], &test[0], &test[1], z,
            acb_mat_entry(tau, 0, 0), prec);
        for (k = 0; k < n2; k++)
        {
            acb_sqr(&test[k], &test[k], prec);
        }

        if (!_acb_vec_overlaps(test, th, n2))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, mat:\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("\n");
            flint_printf("image tau: ");
            acb_printd(acb_mat_entry(tau, 0, 0), 10);
            flint_printf("\n");
            flint_printf("th, test:\n");
            _acb_vec_printd(th, n2, 5);
            flint_printf("\n");
            _acb_vec_printd(test, n2, 5);
        }

        acb_mat_clear(tau);
        fmpz_mat_clear(mat);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, n2);
        _acb_vec_clear(test, n2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
