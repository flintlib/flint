/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_coulomb, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, t, u, eta, l;
        acb_t F1, G1, Hpos1, Hneg1;
        acb_t F2, G2, Hpos2, Hneg2;
        slong prec1, prec2;
        unsigned int mask;

        acb_init(z); acb_init(t); acb_init(u); acb_init(eta); acb_init(l);
        acb_init(F1); acb_init(G1); acb_init(Hpos1); acb_init(Hneg1);
        acb_init(F2); acb_init(G2); acb_init(Hpos2); acb_init(Hneg2);

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        acb_randtest_param(eta, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_randtest_param(l, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_randtest_param(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        acb_randtest_param(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        acb_randtest_param(t, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        acb_add(z, z, t, 1000);
        acb_sub(z, z, t, 1000);

        acb_hypgeom_coulomb(F1, G1, Hpos1, Hneg1, l, eta, z, prec1);
        acb_hypgeom_coulomb(F2, G2, Hpos2, Hneg2, l, eta, z, prec2);

        if (!acb_overlaps(F1, F2) || !acb_overlaps(G1, G2) ||
            !acb_overlaps(Hpos1, Hpos2) || !acb_overlaps(Hneg1, Hneg2))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("l = "); acb_printn(l, 30, 0); flint_printf("\n\n");
            flint_printf("eta = "); acb_printn(eta, 30, 0); flint_printf("\n\n");
            flint_printf("z = "); acb_printn(z, 30, 0); flint_printf("\n\n");
            flint_printf("F1 = "); acb_printn(F1, 30, 0); flint_printf("\n");
            flint_printf("F2 = "); acb_printn(F2, 30, 0); flint_printf("\n\n");
            flint_printf("G1 = "); acb_printn(G1, 30, 0); flint_printf("\n");
            flint_printf("G2 = "); acb_printn(G2, 30, 0); flint_printf("\n\n");
            flint_printf("Hpos1 = "); acb_printn(Hpos1, 30, 0); flint_printf("\n");
            flint_printf("Hpos2 = "); acb_printn(Hpos2, 30, 0); flint_printf("\n\n");
            flint_printf("Hneg1 = "); acb_printn(Hneg1, 30, 0); flint_printf("\n");
            flint_printf("Hneg2 = "); acb_printn(Hneg2, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_mul_onei(t, F1);
        acb_add(t, t, G1, prec1);
        acb_div_onei(u, F1);
        acb_add(u, u, G1, prec1);

        if (!acb_overlaps(t, Hpos1) || !acb_overlaps(u, Hneg1))
        {
            flint_printf("FAIL: complex\n\n");
            flint_printf("l = "); acb_printn(l, 30, 0); flint_printf("\n\n");
            flint_printf("eta = "); acb_printn(eta, 30, 0); flint_printf("\n\n");
            flint_printf("z = "); acb_printn(z, 30, 0); flint_printf("\n\n");
            flint_printf("F1 = "); acb_printn(F1, 30, 0); flint_printf("\n");
            flint_printf("G1 = "); acb_printn(G1, 30, 0); flint_printf("\n");
            flint_printf("Hpos1 = "); acb_printn(Hpos1, 30, 0); flint_printf("\n");
            flint_printf("Hneg1 = "); acb_printn(Hneg1, 30, 0); flint_printf("\n");
            flint_printf("t = "); acb_printn(t, 30, 0); flint_printf("\n");
            flint_printf("u = "); acb_printn(u, 30, 0); flint_printf("\n");
            flint_abort();
        }

        mask = n_randlimb(state);

        if (n_randint(state, 2) == 0)
        {
            acb_randtest_param(t, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
            acb_add(z, z, t, prec2);
            acb_sub(z, z, t, prec2);
        }

        /* also test aliasing */
        acb_set(G2, z);

        acb_hypgeom_coulomb((mask & 1) ? F2 : NULL,
                            (mask & 2) ? G2 : NULL,
                            (mask & 4) ? Hpos2 : NULL,
                            (mask & 8) ? Hneg2 : NULL, l, eta, G2, prec2);

        if (((mask & 1) && !acb_overlaps(F1, F2)) ||
            ((mask & 2) && !acb_overlaps(G1, G2)) ||
            ((mask & 4) && !acb_overlaps(Hpos1, Hpos2)) ||
            ((mask & 8) && !acb_overlaps(Hneg1, Hneg2)))
        {
            flint_printf("FAIL: consistency (mask)\n\n");
            flint_printf("mask = %u\n\n", mask);
            flint_printf("l = "); acb_printn(l, 30, 0); flint_printf("\n\n");
            flint_printf("eta = "); acb_printn(eta, 30, 0); flint_printf("\n\n");
            flint_printf("z = "); acb_printn(z, 30, 0); flint_printf("\n\n");
            flint_printf("F1 = "); acb_printn(F1, 30, 0); flint_printf("\n");
            flint_printf("F2 = "); acb_printn(F2, 30, 0); flint_printf("\n\n");
            flint_printf("G1 = "); acb_printn(G1, 30, 0); flint_printf("\n");
            flint_printf("G2 = "); acb_printn(G2, 30, 0); flint_printf("\n\n");
            flint_printf("Hpos1 = "); acb_printn(Hpos1, 30, 0); flint_printf("\n");
            flint_printf("Hpos2 = "); acb_printn(Hpos2, 30, 0); flint_printf("\n\n");
            flint_printf("Hneg1 = "); acb_printn(Hneg1, 30, 0); flint_printf("\n");
            flint_printf("Hneg2 = "); acb_printn(Hneg2, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        /* Check F_{l-1} G_{l} - F_{l} G_{l-1} = l (l^2 + eta^2)^(1/2) */
        acb_sub_ui(t, l, 1, prec2);
        acb_hypgeom_coulomb(F2, G2, NULL, NULL, t, eta, z, prec2);

        acb_mul_onei(t, eta);
        acb_add(t, t, l, prec2);
        acb_rsqrt(t, t, prec2);
        acb_div_onei(u, eta);
        acb_add(u, u, l, prec2);
        acb_rsqrt(u, u, prec2);
        acb_mul(u, t, u, prec2);
        acb_mul(u, u, l, prec2);

        acb_mul(t, F2, G1, prec2);
        acb_submul(t, F1, G2, prec2);

        if (!acb_overlaps(t, u))
        {
            flint_printf("FAIL: cross-product\n\n");
            flint_printf("l = "); acb_printn(l, 30, 0); flint_printf("\n\n");
            flint_printf("eta = "); acb_printn(eta, 30, 0); flint_printf("\n\n");
            flint_printf("z = "); acb_printn(z, 30, 0); flint_printf("\n\n");
            flint_printf("F1 = "); acb_printn(F1, 30, 0); flint_printf("\n");
            flint_printf("G1 = "); acb_printn(G1, 30, 0); flint_printf("\n");
            flint_printf("F2 = "); acb_printn(F2, 30, 0); flint_printf("\n");
            flint_printf("G2 = "); acb_printn(G2, 30, 0); flint_printf("\n");
            flint_printf("t = "); acb_printn(t, 30, 0); flint_printf("\n");
            flint_printf("u = "); acb_printn(u, 30, 0); flint_printf("\n");
            flint_abort();
        }

        acb_clear(z); acb_clear(t); acb_clear(u); acb_clear(eta); acb_clear(l);
        acb_clear(F1); acb_clear(G1); acb_clear(Hpos1); acb_clear(Hneg1);
        acb_clear(F2); acb_clear(G2); acb_clear(Hpos2); acb_clear(Hneg2);
    }

    TEST_FUNCTION_END(state);
}
