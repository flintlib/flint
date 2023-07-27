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

    flint_printf("transform_sqr_proj....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: inverse matrix gives back the same projective point */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n2 = 1 << (2 * g);
        slong prec = 100;
        slong bits = n_randint(state, 5);
        fmpz_mat_t mat, inv;
        acb_mat_t tau;
        acb_ptr z, th, aux, test;
        acb_t scal;
        slong k;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        fmpz_mat_init(inv, 2 * g, 2 * g);
        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(n2);
        aux = _acb_vec_init(n2);
        test = _acb_vec_init(n2);
        acb_init(scal);

        acb_siegel_randtest_nice(tau, state, prec);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        sp2gz_randtest(mat, state, bits);
        sp2gz_inv(inv, mat);
        acb_theta_ql_all_sqr(test, z, tau, prec);

        acb_theta_transform_sqr_proj(aux, test, mat, prec);
        acb_theta_transform_sqr_proj(th, aux, inv, prec);
        acb_div(scal, &test[0], &th[0], prec);
        _acb_vec_scalar_mul(th, th, n2, scal, prec);

        if (!_acb_vec_overlaps(th, test, n2))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, tau:\n", g);
            acb_mat_printd(tau, 5);
            flint_printf("test, th:\n");
            _acb_vec_printd(test, n2, 5);
            _acb_vec_printd(th, n2, 5);
            flint_abort();
        }

        fmpz_mat_clear(mat);
        fmpz_mat_clear(inv);
        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, n2);
        _acb_vec_clear(aux, n2);
        _acb_vec_clear(test, n2);
        acb_clear(scal);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
