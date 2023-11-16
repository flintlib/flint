/*
    Copyright (C) 2015 Fredrik Johansson

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

TEST_FUNCTION_START(acb_hypgeom_airy, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, t, w;
        acb_t ai1, aip1, bi1, bip1;
        acb_t ai2, aip2, bi2, bip2;
        slong n1, n2, prec1, prec2;
        unsigned int mask;

        acb_init(z); acb_init(t); acb_init(w);
        acb_init(ai1); acb_init(aip1); acb_init(bi1); acb_init(bip1);
        acb_init(ai2); acb_init(aip2); acb_init(bi2); acb_init(bip2);

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        n1 = n_randint(state, 300);
        n2 = n_randint(state, 300);

        acb_randtest_param(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_param(t, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_add(z, z, t, 1000);
        acb_sub(z, z, t, 1000);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_airy_direct(ai1, aip1, bi1, bip1, z, n1, prec1);
                break;
            case 1:
                acb_hypgeom_airy_asymp(ai1, aip1, bi1, bip1, z, n1, prec1);
                break;
            default:
                acb_hypgeom_airy(ai1, aip1, bi1, bip1, z, prec1);
                break;
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_airy_direct(ai2, aip2, bi2, bip2, z, n2, prec2);
                break;
            case 1:
                acb_hypgeom_airy_asymp(ai2, aip2, bi2, bip2, z, n2, prec2);
                break;
            default:
                acb_hypgeom_airy(ai2, aip2, bi2, bip2, z, prec2);
                break;
        }

        if (!acb_overlaps(ai1, ai2))
        {
            flint_printf("FAIL: consistency (Ai)\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("ai1 = "); acb_printd(ai1, 30); flint_printf("\n\n");
            flint_printf("ai2 = "); acb_printd(ai2, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_overlaps(aip1, aip2))
        {
            flint_printf("FAIL: consistency (Ai')\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("aip1 = "); acb_printd(aip1, 30); flint_printf("\n\n");
            flint_printf("aip2 = "); acb_printd(aip2, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_overlaps(bi1, bi2))
        {
            flint_printf("FAIL: consistency (Bi)\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("bi1 = "); acb_printd(bi1, 30); flint_printf("\n\n");
            flint_printf("bi2 = "); acb_printd(bi2, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_overlaps(bip1, bip2))
        {
            flint_printf("FAIL: consistency (Bi')\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("bip1 = "); acb_printd(bip1, 30); flint_printf("\n\n");
            flint_printf("bip2 = "); acb_printd(bip2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_mul(w, ai1, bip1, prec1);
        acb_submul(w, bi1, aip1, prec1);
        acb_const_pi(t, prec1);
        acb_inv(t, t, prec1);

        if (!acb_overlaps(w, t))
        {
            flint_printf("FAIL: wronskian\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("ai1  = "); acb_printd(ai1, 30); flint_printf("\n\n");
            flint_printf("aip1 = "); acb_printd(aip1, 30); flint_printf("\n\n");
            flint_printf("bi1  = "); acb_printd(bi1, 30); flint_printf("\n\n");
            flint_printf("bip1 = "); acb_printd(bip1, 30); flint_printf("\n\n");
            flint_printf("w = ");  acb_printd(w, 30); flint_printf("\n\n");
            flint_abort();
        }

        mask = n_randlimb(state);

        acb_hypgeom_airy((mask & 1) ? ai2 : NULL,
                         (mask & 2) ? aip2 : NULL,
                         (mask & 4) ? bi2 : NULL,
                         (mask & 8) ? bip2 : NULL, z, prec2);

        if (!acb_overlaps(ai1, ai2) || !acb_overlaps(aip1, aip2) ||
            !acb_overlaps(bi1, bi2) || !acb_overlaps(bip1, bip2))
        {
            flint_printf("FAIL: consistency (mask)\n\n");
            flint_printf("mask = %u\n\n", mask);
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("ai1 = "); acb_printd(ai1, 30); flint_printf("\n\n");
            flint_printf("ai2 = "); acb_printd(ai2, 30); flint_printf("\n\n");
            flint_printf("aip1 = "); acb_printd(aip1, 30); flint_printf("\n\n");
            flint_printf("aip2 = "); acb_printd(aip2, 30); flint_printf("\n\n");
            flint_printf("bi1 = "); acb_printd(bi1, 30); flint_printf("\n\n");
            flint_printf("bi2 = "); acb_printd(bi2, 30); flint_printf("\n\n");
            flint_printf("bip1 = "); acb_printd(bip1, 30); flint_printf("\n\n");
            flint_printf("bip2 = "); acb_printd(bip2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z); acb_clear(t); acb_clear(w);
        acb_clear(ai1); acb_clear(aip1); acb_clear(bi1); acb_clear(bip1);
        acb_clear(ai2); acb_clear(aip2); acb_clear(bi2); acb_clear(bip2);
    }

    TEST_FUNCTION_END(state);
}
