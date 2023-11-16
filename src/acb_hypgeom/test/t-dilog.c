/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_dilog, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z0, z1, w0, w1, a;
        slong prec0, prec1;
        int alg0, alg1;

        acb_init(z0);
        acb_init(z1);
        acb_init(w0);
        acb_init(w1);
        acb_init(a);

        prec0 = 2 + n_randint(state, 500);
        prec1 = 2 + n_randint(state, 500);

        acb_randtest(z0, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));

        acb_set(z1, z0);
        if (n_randint(state, 2))
        {
            acb_add(z1, z1, w0, prec0);
            acb_sub(z1, z1, w0, prec0);
        }

        switch ((alg0 = n_randint(state, 15)))
        {
            case 0:
                acb_hypgeom_dilog_zero(w0, z0, prec0);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
                acb_hypgeom_dilog_transform(w0, z0, alg0, prec0);
                break;
            case 8:
                acb_hypgeom_dilog_bernoulli(w0, z0, prec0);
                break;
            case 9:
                acb_hypgeom_dilog_bitburst(w0, a, z0, prec0);
                acb_hypgeom_dilog(a, a, prec0);
                acb_add(w0, w0, a, prec0);
                break;
            default:
                acb_hypgeom_dilog(w0, z0, prec0);
        }

        acb_set(w1, z1);   /* also test aliasing */

        switch ((alg1 = n_randint(state, 15)))
        {
            case 0:
                acb_hypgeom_dilog_zero(w1, w1, prec1);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
                acb_hypgeom_dilog_transform(w1, w1, alg1, prec1);
                break;
            case 8:
                acb_hypgeom_dilog_bernoulli(w1, z1, prec1);
                break;
            case 9:
                acb_hypgeom_dilog_bitburst(w1, a, z1, prec1);
                acb_hypgeom_dilog(a, a, prec1);
                acb_add(w1, w1, a, prec1);
                break;
            default:
                acb_hypgeom_dilog(w1, z1, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("alg0 = %d, alg1 = %d\n\n", alg0, alg1);
            flint_printf("prec0 = %wd, prec1 = %wd\n\n", prec0, prec1);
            flint_printf("z0 = "); acb_printd(z0, 30); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z0);
        acb_clear(z1);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(a);
    }

    TEST_FUNCTION_END(state);
}
