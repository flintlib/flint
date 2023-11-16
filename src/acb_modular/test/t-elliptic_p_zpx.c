/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_elliptic_p_zpx, state)
{
    slong iter;

    /* Test differential equation */
    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau, z;
        acb_ptr g, wp, wp3, wpd, wpd2;
        slong prec, len, i;

        len = 1 + n_randint(state, 15);
        prec = 2 + n_randint(state, 1000);

        acb_init(tau);
        acb_init(z);
        g = _acb_vec_init(2);
        wp = _acb_vec_init(len + 1);
        wp3 = _acb_vec_init(len);
        wpd = _acb_vec_init(len);
        wpd2 = _acb_vec_init(len);

        acb_randtest(tau, state, prec, 10);
        acb_randtest(z, state, prec, 10);

        acb_modular_elliptic_p_zpx(wp, z, tau, len + 1, prec);
        acb_modular_eisenstein(g, tau, 2, prec);
        acb_mul_ui(g, g, 60, prec);
        acb_mul_ui(g + 1, g + 1, 140, prec);

        _acb_poly_derivative(wpd, wp, len + 1, prec);
        _acb_poly_mullow(wpd2, wpd, len, wpd, len, len, prec);
        _acb_poly_pow_ui_trunc_binexp(wp3, wp, len, 3, len, prec);
        _acb_vec_scalar_mul_ui(wp3, wp3, len, 4, prec);
        _acb_vec_scalar_submul(wp3, wp, len, g, prec);
        acb_sub(wp3, wp3, g + 1, prec);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(wpd2 + i, wp3 + i))
            {
                flint_printf("FAIL (overlap)\n");
                flint_printf("i = %wd  len = %wd  prec = %wd\n\n", i, len, prec);
                flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
                flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
                flint_printf("wp = "); acb_printd(wp + i, 15); flint_printf("\n\n");
                flint_printf("wpd = "); acb_printd(wpd + i, 15); flint_printf("\n\n");
                flint_printf("wp3 = "); acb_printd(wp3 + i, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(tau);
        acb_clear(z);
        _acb_vec_clear(g, 2);
        _acb_vec_clear(wp, len + 1);
        _acb_vec_clear(wp3, len);
        _acb_vec_clear(wpd, len);
        _acb_vec_clear(wpd2, len);
    }

    /* Consistency test */
    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau, z;
        acb_ptr wp1, wp2;
        slong prec1, prec2, len1, len2, i;

        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_init(tau);
        acb_init(z);
        wp1 = _acb_vec_init(len1);
        wp2 = _acb_vec_init(len2);

        acb_randtest(tau, state, prec1, 10);
        acb_randtest(z, state, prec1, 10);

        acb_modular_elliptic_p_zpx(wp1, z, tau, len1, prec1);
        acb_modular_elliptic_p_zpx(wp2, z, tau, len2, prec2);

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            if (!acb_overlaps(wp1 + i, wp2 + i))
            {
                flint_printf("FAIL (overlap)\n");
                flint_printf("i = %wd len1 = %wd len2 = %wd\n\n", i, len1, len2);
                flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
                flint_printf("wp1 = "); acb_printd(wp1 + i, 15); flint_printf("\n\n");
                flint_printf("wp2 = "); acb_printd(wp2 + i, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(tau);
        acb_clear(z);
        _acb_vec_clear(wp1, len1);
        _acb_vec_clear(wp2, len2);
    }

    TEST_FUNCTION_END(state);
}
