/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_sextic_chi5, state)
{
    slong iter;

    /* Test: discriminant of sextic is chi10 = chi5^2 */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 1000);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, th2, roots, mf;
        acb_poly_t f;
        acb_t chi5, d, t;
        slong nb, k, j;
        int res = 1;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        roots = _acb_vec_init(6);
        mf = _acb_vec_init(4);
        acb_poly_init(f);
        acb_init(chi5);
        acb_init(d);
        acb_init(t);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);

        acb_theta_g2_sextic_chi5(f, chi5, tau, prec);

        for (k = 0; k <= 6; k++)
        {
            acb_poly_get_coeff_acb(t, f, k);
            res = res && acb_is_finite(t);
        }

        if ((!res && !acb_contains_zero(chi5))
            || !acb_is_finite(chi5))
        {
            flint_printf("FAIL (finite)\n");
            acb_poly_printd(f, 5);
            flint_printf("\n");
            flint_abort();
        }

        nb = acb_poly_find_roots(roots, f, NULL, 0, prec);

        if (nb == 6)
        {
            acb_one(d);
            for (k = 0; k < 6; k++)
            {
                for (j = k + 1; j < 6; j++)
                {
                    acb_sub(t, &roots[k], &roots[j], prec);
                    acb_mul(d, d, t, prec);
                }
            }
            acb_sqr(d, d, prec);
            acb_poly_get_coeff_acb(t, f, 6);
            acb_pow_ui(t, t, 10, prec);
            acb_mul(d, d, t, prec);
            acb_mul_2exp_si(d, d, -12);

            acb_theta_all(th2, z, tau, 1, prec);
            acb_theta_g2_even_weight(&mf[0], &mf[1], &mf[2], &mf[3], th2, prec);
            acb_sqr(chi5, chi5, prec);

            if (!acb_overlaps(d, &mf[2])
                || !acb_overlaps(chi5, &mf[2])
                || !_acb_vec_is_finite(mf, 4)
                || !acb_is_finite(d))
            {
                flint_printf("FAIL\n");
                flint_printf("roots, discr, chi10, chi5^2:\n");
                _acb_vec_printd(roots, 6, 5);
                acb_printd(d, 5);
                flint_printf("\n");
                acb_printd(&mf[2], 5);
                flint_printf("\n");
                acb_printd(chi5, 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        _acb_vec_clear(roots, 6);
        _acb_vec_clear(mf, 4);
        acb_poly_clear(f);
        acb_clear(chi5);
        acb_clear(d);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
