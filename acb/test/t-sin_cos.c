/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("sin_cos....");
    fflush(stdout);

    flint_randinit(state);

    /* check sin(a+b) = cos(b)*sin(a) + cos(a)*sin(b) */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, cosa, sina, cosb, sinb;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(cosa);
        acb_init(sina);
        acb_init(cosb);
        acb_init(sinb);

        acb_randtest(a, state, 1 + n_randint(state, 200), 3);
        acb_randtest(b, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        acb_add(c, a, b, prec);
        acb_sin(c, c, prec);

        acb_sin_cos(sina, cosa, a, prec);
        acb_sin_cos(sinb, cosb, b, prec);
        acb_mul(cosb, cosb, sina, prec);
        acb_mul(cosa, cosa, sinb, prec);
        acb_add(d, cosa, cosb, prec);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: sin(a+b) = cos(b)*sin(a) + cos(a)*sin(b)\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(cosa);
        acb_clear(sina);
        acb_clear(cosb);
        acb_clear(sinb);
    }

    /* check cos(a+b) = cos(b)*cos(a) - sin(a)*sin(b) */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, cosa, sina, cosb, sinb;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(cosa);
        acb_init(sina);
        acb_init(cosb);
        acb_init(sinb);

        acb_randtest(a, state, 1 + n_randint(state, 200), 3);
        acb_randtest(b, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        acb_add(c, a, b, prec);
        acb_cos(c, c, prec);

        acb_sin_cos(sina, cosa, a, prec);
        acb_sin_cos(sinb, cosb, b, prec);
        acb_mul(cosa, cosa, cosb, prec);
        acb_mul(sina, sina, sinb, prec);
        acb_sub(d, cosa, sina, prec);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: cos(a+b) = cos(b)*cos(a) - sin(a)*sin(b)\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(cosa);
        acb_clear(sina);
        acb_clear(cosb);
        acb_clear(sinb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
