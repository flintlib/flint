/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"
#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("fast_addmul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y, z, z2, w;
        mag_t xb, yb, zb;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(z2);
        fmpr_init(w);

        mag_init(xb);
        mag_init(yb);
        mag_init(zb);

        mag_randtest(xb, state, 15);
        mag_randtest(yb, state, 15);
        mag_randtest(zb, state, 15);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);
        mag_get_fmpr(z, zb);

        fmpr_addmul(z, x, y, MAG_BITS + 10, FMPR_RND_DOWN);
        fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(z2, z2, -10);

        mag_fast_addmul(zb, xb, yb);
        mag_get_fmpr(w, zb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)
        MAG_CHECK_BITS(zb)

        if (!(fmpr_cmpabs(z, w) <= 0 && fmpr_cmpabs(w, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); fmpr_printd(z, 15); flint_printf("\n\n");
            flint_printf("w = "); fmpr_printd(w, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_set(z, x);
        fmpr_addmul(z, z, y, MAG_BITS + 10, FMPR_RND_DOWN);
        fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(z2, z2, -10);

        mag_set(zb, xb);
        mag_fast_addmul(zb, zb, yb);
        mag_get_fmpr(w, zb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)
        MAG_CHECK_BITS(zb)

        if (!(fmpr_cmpabs(z, w) <= 0 && fmpr_cmpabs(w, z2) <= 0))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_printf("x = "); fmpr_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); fmpr_printd(z, 15); flint_printf("\n\n");
            flint_printf("w = "); fmpr_printd(w, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_set(z, y);
        fmpr_addmul(z, x, z, MAG_BITS + 10, FMPR_RND_DOWN);
        fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(z2, z2, -10);

        mag_set(zb, yb);
        mag_fast_addmul(zb, xb, zb);
        mag_get_fmpr(w, zb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)
        MAG_CHECK_BITS(zb)

        if (!(fmpr_cmpabs(z, w) <= 0 && fmpr_cmpabs(w, z2) <= 0))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
            flint_printf("x = "); fmpr_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); fmpr_printd(z, 15); flint_printf("\n\n");
            flint_printf("w = "); fmpr_printd(w, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_set(z, x);
        fmpr_addmul(z, z, z, MAG_BITS + 10, FMPR_RND_DOWN);
        fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(z2, z2, -10);

        mag_set(zb, xb);
        mag_fast_addmul(zb, zb, zb);
        mag_get_fmpr(w, zb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)
        MAG_CHECK_BITS(zb)

        if (!(fmpr_cmpabs(z, w) <= 0 && fmpr_cmpabs(w, z2) <= 0))
        {
            flint_printf("FAIL (aliasing 3)\n\n");
            flint_printf("x = "); fmpr_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); fmpr_printd(z, 15); flint_printf("\n\n");
            flint_printf("w = "); fmpr_printd(w, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(z2);
        fmpr_clear(w);

        mag_clear(xb);
        mag_clear(yb);
        mag_clear(zb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

