/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_fresnel, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, z2, s, c, u, v;
        slong prec1, prec2;
        int normalized;

        prec1 = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        acb_init(z);
        acb_init(z2);
        acb_init(s);
        acb_init(c);
        acb_init(u);
        acb_init(v);

        acb_randtest_special(z, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest_special(s, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest_special(c, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));

        normalized = n_randint(state, 2);

        /* test S(z) + i C(z) = sqrt(pi/2) (1+i)/2 erf((1+i)/sqrt(2) z) */

        /* u = rhs */
        acb_onei(u);
        acb_sqrt(u, u, prec1);
        acb_mul(u, u, z, prec1);
        acb_hypgeom_erf(u, u, prec1);
        acb_mul_onei(v, u);
        acb_add(u, u, v, prec1);
        acb_mul_2exp_si(u, u, -1);
        acb_const_pi(v, prec1);
        acb_mul_2exp_si(v, v, -1);
        acb_sqrt(v, v, prec1);
        acb_mul(u, u, v, prec1);

        if (normalized)
        {
            acb_const_pi(v, prec2);
            acb_mul_2exp_si(v, v, -1);
            acb_sqrt(v, v, prec2);
            acb_div(z2, z, v, prec2);
        }
        else
        {
            acb_set(z2, z);
        }

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_fresnel(s, c, z2, normalized, prec2);
                break;
            case 1:
                acb_hypgeom_fresnel(s, NULL, z2, normalized, prec2);
                acb_hypgeom_fresnel(NULL, c, z2, normalized, prec2);
                break;
            case 2:
                acb_set(s, z2);
                acb_hypgeom_fresnel(s, c, s, normalized, prec2);
                break;
            case 3:
                acb_set(c, z2);
                acb_hypgeom_fresnel(s, c, c, normalized, prec2);
                break;
            default:
                acb_hypgeom_fresnel(s, c, z2, normalized, prec2);
        }

        if (normalized)
        {
            acb_mul(s, s, v, prec2);
            acb_mul(c, c, v, prec2);
        }

        acb_mul_onei(v, c);
        acb_add(v, v, s, prec2);

        if (!acb_overlaps(u, v))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_printf("u = "); acb_printd(u, 30); flint_printf("\n\n");
            flint_printf("v = "); acb_printd(v, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(z2);
        acb_clear(s);
        acb_clear(c);
        acb_clear(u);
        acb_clear(v);
    }

    TEST_FUNCTION_END(state);
}
